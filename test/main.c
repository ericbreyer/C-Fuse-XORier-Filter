//------------------------------------------------------------------------------
// Information
//------------------------------------------------------------------------------

/**
 * @file main.c
 * @author Eric Breyer (eric.breyer@gmail.com) [https://github.com/ericbreyer]
 * @brief Main file for testing and profiling fuse XORier lookup table
 * @version 1.0
 * @date 2023-12-13
 *
 * @copyright Copyright (c) 2023, released under GNU Public Licence 3.0 or later
 *
 */
// SPDX-License-Identifier: GNU-3.0-or-later

//------------------------------------------------------------------------------
// Includes
//------------------------------------------------------------------------------

#include <math.h>
#include <omp.h>
#include <stdint.h> // printing
#include <stdio.h>  // printing
#include <stdlib.h> // malloc
#include <sys/resource.h>
#include <sys/stat.h> /* For mode constants */
#include <time.h>

#include <fuse_XORier_lookup_table.h>
#include "gnuplot_i.h"
#include "hash_table.h"
#include "linked_hash_table.h"
#include "rbtree.h"

//------------------------------------------------------------------------------
// Defines
//------------------------------------------------------------------------------

/**
 * @brief The maximum number of threads to use for parallelization
 * 
 */
#define MAX_THREADS 4

//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------

void generate_hash_advantage_plot(char * outfile, int start, int end, int step, double c,
                                  int k) {
    const int NPOINTS = (end - start) / step + 1;
    struct {
        double time;
        size_t mem;
        size_t peak_mem;
        size_t kvmem;
        size_t n
    } *testData = malloc(sizeof *testData * NPOINTS * 2);
    size_t ns[NPOINTS * 2];
    int i = 0;
    for (int n = start; n <= end; n += step, ++i) {
        testData[i].n = ns[i] = n;
        testData[i + NPOINTS].n = ns[i + NPOINTS] = n;
    }

#pragma omp parallel for num_threads(MAX_THREADS) schedule(static, 1)
    for (i = 0; i < NPOINTS * 2; ++i) {

        int n = ns[i];
        fprintf(stderr, "Start Fuse Bloomier Filter with n=%d\n", n);
        SizedPointer *keys = malloc(sizeof *keys * n);
        void **vals = malloc(sizeof *vals * n);
        for (int i = 0; i < n; ++i) {
            char *intStr = malloc(sizeof *intStr * (i % 10 + 2 + 15));
            if (intStr == NULL) {
                perror("malloc");
                exit(1);
            }
            int len = sprintf(intStr, "%d", i);
            keys[i].ptr = intStr;
            keys[i].size = len;
            vals[i] = (void *)(uintptr_t)((i & 1) + 1);
        }

        size_t kvmem = n * sizeof(void *);
        for (int j = 0; j < n; ++j) {
            kvmem += keys[j].size * sizeof(char);
        }
        testData[i].kvmem = kvmem;

        // struct BloomierFilter *fxlt = construct_bloomierFilter(k, v, n, 2.2,
        // 2,  0); struct BloomierFilter *fxlt = construct_bloomierFilter(k,
        // v, n, 1.23, 3,  0); struct BloomierFilter *fxlt =
        // construct_bloomierFilter(k, v, n, 1.15, 3,  nfor3); struct
        // BloomierFilter *fxlt = construct_bloomierFilter(k, v, n, 1.075, 4,
        // nfor4);

        struct timespec start, finish;

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct fuseXORierLookupTable *fxlt = build_fuseXORierLT(
            n, keys, vals, c, k, i >= NPOINTS ? FXLT_FLAG_CACHE_HASHES : 0);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);

        testData[i].time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;

        struct footprint fp2 = fuseXORierLT_getSize(fxlt);
        testData[i].mem = fp2.size;
        testData[i].peak_mem = fp2.peak_size;

        delete_fuseXORierLT(fxlt);
        fprintf(stderr, "Fuse Bloomier Filter with n=%d c=%f done in %fs\n", n,
                c, testData[i].time);
        for (int i = 0; i < n; ++i) {
            free(keys[i].ptr);
        }
        free(keys);
        free(vals);
    }

    gnuplot_ctrl *h1 = gnuplot_init();
    // gnuplot_setstyle(h1, "lines");
    gnuplot_cmd(h1, "set terminal png size 1000,500\n");
    char * setoutputcmd = malloc(sizeof *setoutputcmd * (strlen(outfile) + 20));
    sprintf(setoutputcmd, "set output \"%s\"\n", outfile);
    gnuplot_cmd(h1, setoutputcmd);
    gnuplot_cmd(h1, "set multiplot layout 1,2\n");

    // print the n,time pairs to temp.dat
    FILE *temp = fopen("nocache.dat", "w");
    for (int i = 0; i < NPOINTS; ++i) {
        fprintf(temp, "%f %f %f %f %f\n", (double)testData[i].n,
                testData[i].time, (double)testData[i].mem / 1000000,
                (double)testData[i].peak_mem / 1000000,
                (double)testData[i].kvmem / 1000000);
    }
    fflush(temp);
    fclose(temp);
    temp = fopen("cache.dat", "w");
    for (int i = NPOINTS; i < NPOINTS * 2; ++i) {
        fprintf(temp, "%f %f %f %f\n", (double)testData[i].n, testData[i].time,
                (double)testData[i].mem / 1000000,
                (double)testData[i].peak_mem / 1000000);
    }
    fflush(temp);
    fclose(temp);

    gnuplot_cmd(h1, "set xlabel \"number of keys\"\n");
    gnuplot_cmd(h1, "set ylabel \"build time\"\n");
    gnuplot_cmd(h1, "set title \"Build Time vs. Number of Keys\"\n");
    gnuplot_cmd(h1, "set colorbox vertical origin screen 0.9, 0.2 size screen "
                    "0.05, 0.6 front  noinvert bdefault\n");
    gnuplot_cmd(h1, "plot \"nocache.dat\" using 1:2 linecolor rgb '#dd181f' "
                    "linetype 1 linewidth 2 pt 7 ps 1 title 'time',"
                    " \"cache.dat\"   using 1:2 linecolor rgb '#2CB916' "
                    "linetype 1 linewidth 2 pt 7 ps 1 title 'cache time'\n");
    gnuplot_cmd(h1, "set xlabel \"number of keys\"\n");
    gnuplot_cmd(h1, "set ylabel \"storage (mb)\"\n");
    gnuplot_cmd(h1, "set title \"Space vs. Number of Keys\"\n");
    gnuplot_cmd(h1, "set colorbox vertical origin screen 0.9, 0.2 size screen "
                    "0.05, 0.6 front  noinvert bdefault\n");
    gnuplot_cmd(h1,
                "plot \"nocache.dat\" using 1:3 linecolor rgb '#A5B2BD' "
                "linetype 1 linewidth 2 pt 7 ps 1 title 'memory',"
                "\"nocache.dat\" using 1:5 linecolor rgb '#6B737A' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'raw k/v memory',"
                "\"nocache.dat\" using 1:4 linecolor rgb '#dd181f' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'peak memory',"
                "\"cache.dat\"   using 1:4 linecolor rgb '#2CB916' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'cache peak memory'\n");

    gnuplot_close(h1);
}
void generate_FXLT_vs_other_maps_plot(char * outfile, int start, int end, int step, double c,
                                      int k) {
    struct timeMem {
        double build_time;
        double lookup_time;
        size_t mem;
    };
    struct timeMemComp {
        struct timeMem CFXLT;
        struct timeMem FXLT;
        struct timeMem HT;
        struct timeMem LHT;
        // struct timeMem RBT;
    };

    const int NPOINTS = abs(end - start) / step + 1;
    struct timeMemComp *results = malloc(sizeof *results * NPOINTS);

#pragma omp parallel for num_threads(MAX_THREADS) schedule(static, 1)
    for (int i = 0; i < NPOINTS; ++i) {
        int n = start + i * step;
        fprintf(stderr, "Start Test with n=%d on thread %d\n", n,
                omp_get_thread_num());

        SizedPointer *keys = malloc(sizeof *keys * n);
        void **vals = malloc(sizeof *vals * n);
        for (int i = 0; i < n; ++i) {
            char *intStr = malloc(sizeof *intStr * (i % 10 + 10));
            int len = sprintf(intStr, "%d long", i);
            keys[i].ptr = intStr;
            keys[i].size = len;
            vals[i] = (void *)(uintptr_t)((i & 1) + 1);
        }

        size_t kvmem = n * sizeof(void *);
        for (int j = 0; j < n; ++j) {
            kvmem += keys[j].size * sizeof(char);
        }

        struct timespec start, finish;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct fuseXORierLookupTable *fxlt =
            build_fuseXORierLT(n, keys, vals, c, k, FXLT_FLAG_CACHE_HASHES);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].CFXLT.build_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        for (int i = 0; i < n; ++i) {
            fuseXORierLT_lookup(fxlt, keys[i]);
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);

        results[i].CFXLT.lookup_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;

        results[i].CFXLT.mem = fuseXORierLT_getSize(fxlt).size;

        // same thing but no chaching
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct fuseXORierLookupTable *bf2 =
            build_fuseXORierLT(n, keys, vals, c, k, 0);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].FXLT.build_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        for (int i = 0; i < n; ++i) {
            fuseXORierLT_lookup(bf2, keys[i]);
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].FXLT.lookup_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        results[i].FXLT.mem = fuseXORierLT_getSize(bf2).size;

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);

        struct HashTable *t = construct_HashTable();
        for (int i = 0; i < n; ++i) {
            HashTable_insert(t, keys[i], vals[i]);
        }

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].HT.build_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        for (int i = 0; i < n; ++i) {
            HashTable_getValue(t, keys[i]);
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].HT.lookup_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        results[i].HT.mem = HashTable_getSize(t) + kvmem;

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct LinkedHashTable *lt = construct_LinkedHashTable();

        for (int i = 0; i < n; ++i) {
            LinkedHashTable_insert(lt, keys[i].ptr, vals[i]);
        }

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].LHT.build_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        for (int i = 0; i < n; ++i) {
            LinkedHashTable_getValue(lt, keys[i].ptr);
        }
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        results[i].LHT.lookup_time =
            (finish.tv_sec - start.tv_sec) +
            (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        results[i].LHT.mem = LinkedHashTable_getSize(lt) + kvmem;

        // clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        // struct rbtree *rbt = construct_rbtree();

        // for (int i = 0; i < n/ 100; ++i) {
        //     rbtree_insert(rbt, keys[i].ptr, vals[i]);
        // }
        // clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        // results[i].RBT.build_time =
        //     (finish.tv_sec - start.tv_sec) +
        //     (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;

        // clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        // for (int i = 0; i < n/ 100; ++i) {
        //     rbtree_getValue(rbt, keys[i].ptr);
        // }
        // clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        // results[i].RBT.lookup_time =
        //     (finish.tv_sec - start.tv_sec) +
        //     (double)(finish.tv_nsec - start.tv_nsec) * 100 / 1000000000.0;

        // results[i].RBT.mem = (rbtree_getSize(rbt)* 100) + kvmem;

        delete_fuseXORierLT(fxlt);
        delete_fuseXORierLT(bf2);
        destroy_HashTable(t);
        destroy_LinkedHashTable(lt);
        // destroy_rbtree(rbt);
        for (int i = 0; i < n; ++i) {
            free(keys[i].ptr);
        }
        free(keys);
        free(vals);
    }
    fprintf(stderr, "Generating Plot\n");

    gnuplot_ctrl *h1 = gnuplot_init();
    gnuplot_cmd(h1, "set terminal png size 1000,1000\n");
    char * setoutputcmd = malloc(sizeof *setoutputcmd * (strlen(outfile) + 20));
    sprintf(setoutputcmd, "set output \"%s\"\n", outfile);
    gnuplot_cmd(h1, setoutputcmd);
    gnuplot_cmd(h1, "set multiplot layout 2,2\n");
    gnuplot_cmd(h1, "set key at screen 0.975, screen 0.350 font ',16' box ls 1 "
                    "lw 2 lc 'black' spacing 2\n");

    // print the n,time pairs to temp.dat
    FILE *temp = fopen("temp.dat", "w");
    for (int i = 0; i < NPOINTS; ++i) {
        // print all these to file
        int n = start + i * step;
        fprintf(temp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
                (double)n, (double)results[i].CFXLT.mem / 1000000,
                (double)results[i].CFXLT.build_time,
                (double)results[i].CFXLT.lookup_time * 1000000 / n,
                (double)results[i].FXLT.mem / 1000000,
                (double)results[i].FXLT.build_time,
                (double)results[i].FXLT.lookup_time * 1000000 / n,
                (double)results[i].HT.mem / 1000000,
                (double)results[i].HT.build_time,
                (double)results[i].HT.lookup_time * 1000000 / n,
                (double)results[i].LHT.mem / 1000000,
                (double)results[i].LHT.build_time,
                (double)results[i].LHT.lookup_time * 1000000 / n,
                // (double)results[i].RBT.mem / 1000000,
                // (double)results[i].RBT.build_time,
                // (double)results[i].RBT.lookup_time * 1000000 / n
                0.0, 0.0, 0.0);
    }
    fflush(temp);
    fclose(temp);

    gnuplot_cmd(h1, "set xlabel \"number of keys\"\n");
    gnuplot_cmd(h1, "set ylabel \"build time (s)\"\n");
    gnuplot_cmd(h1, "set title \"Build Time (s) vs. Number of Keys\"\n");
    // gnuplot_cmd(h1, "set colorbox vertical origin screen 0.9, 0.2 size screen
    // "
    //                 "0.05, 0.6 front  noinvert bdefault\n");
    // gnuplot_cmd(h1, "set nokey\n");
    gnuplot_cmd(h1,
                "plot \"temp.dat\" using 1:3 linecolor rgb '#dd181f' "
                "linetype 1 linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup "
                "Table w/ Cache',"
                "\"temp.dat\" using 1:6 linecolor rgb '#2CB916' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup Table',"
                "\"temp.dat\" using 1:9 linecolor rgb '#6B737A' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Hash Table (Open Addressing)',"
                "\"temp.dat\" using 1:12 linecolor rgb '#A5B2BD' linetype "
                "1 linewidth 2 pt 7 ps 1 title 'Hash Table (Chaining)'"
                // "\"temp.dat\" using 1:15 linecolor rgb '#000000' linetype "
                // "1 linewidth 2 pt 7 ps 1 title 'Red Black Tree'\n"
                "\n");

    gnuplot_cmd(h1, "set xlabel \"number of keys\"\n");
    gnuplot_cmd(h1, "set ylabel \"average lookup time (µs)\"\n");
    gnuplot_cmd(h1,
                "set title \"Average Lookup Time (µs) vs. Number of Keys\"\n");
    gnuplot_cmd(h1, "set colorbox vertical origin screen 0.9, 0.2 size screen "
                    "0.05, 0.6 front  noinvert bdefault\n");
    gnuplot_cmd(h1,
                "plot \"temp.dat\" using 1:4 linecolor rgb '#dd181f' "
                "linetype 1 linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup "
                "Table w/ Cache',"
                "\"temp.dat\" using 1:7 linecolor rgb '#2CB916' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup Table',"
                "\"temp.dat\" using 1:10 linecolor rgb '#6B737A' linetype "
                "1 linewidth 2 pt 7 ps 1 title 'Hash Table (Open Addressing)',"
                "\"temp.dat\" using 1:13 linecolor rgb '#A5B2BD' linetype "
                "1 linewidth 2 pt 7 ps 1 title 'Hash Table (Chaining)'"
                // "\"temp.dat\" using 1:16 linecolor rgb '#000000' linetype "
                // "1 linewidth 2 pt 7 ps 1 title 'Red Black Tree'\n"
                "\n");

    gnuplot_cmd(h1, "set xlabel \"number of keys (mb)\"\n");
    gnuplot_cmd(h1, "set ylabel \"memory (mb)\"\n");
    gnuplot_cmd(h1, "set title \"Memory vs. Number of Keys\"\n");
    gnuplot_cmd(h1, "set colorbox vertical origin screen 0.9, 0.2 size screen "
                    "0.05, 0.6 front  noinvert bdefault\n");
    gnuplot_cmd(h1,
                "plot \"temp.dat\" using 1:2 linecolor rgb '#dd181f' "
                "linetype 1 linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup "
                "Table w/ Cache',"
                "\"temp.dat\" using 1:5 linecolor rgb '#2CB916' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup Table',"
                "\"temp.dat\" using 1:8 linecolor rgb '#6B737A' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Hash Table (Open Addressing)',"
                "\"temp.dat\" using 1:11 linecolor rgb '#A5B2BD' linetype "
                "1 linewidth 2 pt 7 ps 1 title 'Hash Table (Chaining)'"
                // "\"temp.dat\" using 1:14 linecolor rgb '#000000' linetype "
                // "1 linewidth 2 pt 7 ps 1 title 'Red Black Tree'\n"
                "\n");

    // gnuplot_plot_xy(h1, ns, fps, NPOINTS, "fp rate");
    gnuplot_close(h1);
}

// compare a fuse bloomier filter to one without windowing
void generate_spacial_coupling_advantage_plot(char * outfile, double start, double end,
                                              double step, double n) {
    struct timeComp {
        double FXLT3;
        double FXLT4;
        double BF;
    };

    const int NPOINTS = fabs(end - start) / step + 1;
    struct timeComp *results = malloc(sizeof *results * NPOINTS);
#pragma omp parallel for num_threads(MAX_THREADS) schedule(static, 1)
    for (int i = 0; i < NPOINTS; ++i) {
        double c = start + i * step;
        fprintf(stderr, "Start Test with c=%f on thread %d\n", c,
                omp_get_thread_num());
        SizedPointer *keys = malloc(sizeof *keys * n);
        void **vals = malloc(sizeof *vals * n);
        for (int i = 0; i < n; ++i) {
            char *intStr = malloc(sizeof *intStr * (i % 10 + 10));
            int len = sprintf(intStr, "%d long", i);
            keys[i].ptr = intStr;
            keys[i].size = len;
            vals[i] = (void *)(uintptr_t)((i & 1) + 1);
        }

        size_t kvmem = n * sizeof(void *);
        for (int j = 0; j < n; ++j) {
            kvmem += keys[j].size * sizeof(char);
        }
        int k = 3;

        struct timespec start, finish;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct fuseXORierLookupTable *fxlt =
            build_fuseXORierLT(n, keys, vals, c, k, 0);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        if (fxlt != NULL) {
            results[i].FXLT3 =
                (finish.tv_sec - start.tv_sec) +
                (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            delete_fuseXORierLT(fxlt);
            printf("FXLT3 with c=%f done in %fs\n", c, results[i].FXLT3);

        } else {
            results[i].FXLT3 = 0.1;
            printf("FXLT3 with c=%f failed\n", c);
        }

        k = 4;

        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct fuseXORierLookupTable *bf3 =
            build_fuseXORierLT(n, keys, vals, c, k, 0);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        if (bf3 != NULL) {
            results[i].FXLT4 =
                (finish.tv_sec - start.tv_sec) +
                (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            delete_fuseXORierLT(bf3);
            printf("FXLT4 with c=%f done in %fs\n", c, results[i].FXLT4);

        } else {
            results[i].FXLT4 = 0.05;
            printf("FXLT4 with c=%f failed\n", c);
        }

        k = 3;
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
        struct fuseXORierLookupTable *bf2 = build_fuseXORierLT(
            n, keys, vals, c, k, FXLT_FLAG_NO_SPATIAL_COUPLING);
        clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
        if (bf2 != NULL) {
            results[i].BF =
                (finish.tv_sec - start.tv_sec) +
                (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0;
            delete_fuseXORierLT(bf2);
            printf("BF with c=%f done in %fs\n", c, results[i].BF);

        } else {
            results[i].BF = 0;
            printf("BF with c=%f failed\n", c);
        }

        for (int i = 0; i < n; ++i) {
            free(keys[i].ptr);
        }
        free(keys);
        free(vals);
    }
    fprintf(stderr, "Generating Plot\n");
    gnuplot_ctrl *h1 = gnuplot_init();
    gnuplot_cmd(h1, "set terminal png size 500,500\n");
        char * setoutputcmd = malloc(sizeof *setoutputcmd * (strlen(outfile) + 20));
    sprintf(setoutputcmd, "set output \"%s\"\n", outfile);
    gnuplot_cmd(h1, setoutputcmd);
    // gnuplot_cmd(h1, "set multiplot layout 2,2\n");
    gnuplot_cmd(h1, "set key at screen 0.975, screen 0.350 font ',8' box ls 1 "
                    "lw 2 lc 'black' spacing 2\n");

    // print the n,time pairs to temp.dat
    FILE *temp = fopen("temp.dat", "w");
    for (int i = 0; i < NPOINTS; ++i) {
        // print all these to file
        double c = start + i * step;
        fprintf(temp, "%f %f %f %f\n", (double)c, (double)results[i].FXLT3,
                (double)results[i].FXLT4, (double)results[i].BF);
    }
    fflush(temp);
    fclose(temp);

    gnuplot_cmd(h1, "set xlabel \"c\"\n");
    gnuplot_cmd(h1, "set ylabel \"build time (s)\"\n");
    gnuplot_cmd(h1, "set title \"Build Time (s) vs. c (where m = c * n, n = "
                    "1,000,000)\"\n");
    gnuplot_cmd(h1,
                "plot \"temp.dat\" using 1:2 linecolor rgb '#dd181f' "
                "linetype 1 linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup "
                "Table k = 3',"
                "\"temp.dat\" using 1:3 linecolor rgb '#2CB916' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Fuse XORier Lookup Table k = 4',"
                "\"temp.dat\" using 1:4 linecolor rgb '#6B737A' linetype 1 "
                "linewidth 2 pt 7 ps 1 title 'Bloomier Filter'\n");
    gnuplot_close(h1);
    free(results);
}

static double thread_progress[MAX_THREADS];
void flush_thread_progress() {
    fprintf(stderr, "\r");
    for (int i = 0; i < MAX_THREADS; ++i) {
        fprintf(stderr, "%d: %.2f%%, ", i, thread_progress[i] * 100);
    }
    fflush(stdout);
}
void print_thread_progress(int thread_num, int n, int NPOINTS) {
    thread_progress[thread_num] = (double)n / NPOINTS;
    flush_thread_progress();
}
void reset_thread_progress() {
    for (int i = 0; i < MAX_THREADS; ++i) {
        thread_progress[i] = 0.0;
    }
}

bool **generate_FXLT_build_success_plot_helper(int nStart, int nEnd, int nStep,
                                               double cStart, double cEnd,
                                               double cStep, int k) {
    const int NPOINTS = (nEnd - nStart) / nStep + 1;
    const int CPOINTS = (cEnd - cStart) / cStep + 1;

    bool(*built)[CPOINTS] = malloc(sizeof *built * NPOINTS);

    for (int i = 0; i < NPOINTS; ++i) {
        for (int j = 0; j < CPOINTS; ++j) {
            built[i][j] = 0.0;
        }
    }
    reset_thread_progress();
#pragma omp parallel for ordered num_threads(MAX_THREADS) schedule(static, 1)
    for (int nIdx = 0; nIdx < NPOINTS; ++nIdx) {
        int n = nStart + nIdx * nStep;
        SizedPointer *keys = malloc(sizeof *keys * n);
        void **vals = malloc(sizeof *vals * n);
        for (int i = 0; i < n; ++i) {
            char *intStr = malloc(sizeof *intStr * (i % 10 + 10));
            int len = sprintf(intStr, "%d long", i);
            keys[i].ptr = intStr;
            keys[i].size = len;
            vals[i] = (void *)(uintptr_t)((i & 1) + 1);
        }

        size_t kvmem = n * sizeof(void *);
        for (int j = 0; j < n; ++j) {
            kvmem += keys[j].size * sizeof(char);
        }

        for (int cIdx = 0; cIdx < CPOINTS; ++cIdx) {
            // fprintf(stderr, "Start Test with c=%f n=%d on thread %d\n",
            //         cStart + cIdx * cStep, n, omp_get_thread_num());
            double c = cStart + cIdx * cStep;
            struct timespec start, finish;
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
            struct fuseXORierLookupTable *bf2 =
                build_fuseXORierLT(n, keys, vals, c, k, 0);
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
            if (bf2 != NULL) {
                built[nIdx][cIdx] = true;
                delete_fuseXORierLT(bf2);
                // printf("BF with c=%f n=%d done in %fs\n", c, n,
                // built[nIdx][cIdx]);

            } else {
                built[nIdx][cIdx] = false;
                // printf("BF with c=%f  n = %d, failed\n", c, n);
            }
            // print done with c
            // printf("\rDone with c = %.03f      ", c);
            print_thread_progress(omp_get_thread_num(), cIdx + 1, CPOINTS);
            fflush(stdout);
        }
        for (int i = 0; i < n; ++i) {
            free(keys[i].ptr);
        }
        free(keys);
        free(vals);
        // print done with n
        // #pragma omp ordered
        printf("\rDone with n = %d in thread %d                                "
               " \n",
               n, omp_get_thread_num());
        flush_thread_progress();
    }
    return built;
}

void generate_FXLT_build_success_plot(char * outfile, double cStart, double cEnd, double cStep,
                                      int k) {
#define numIntervals 3
    int nStarts[numIntervals] = {100000, 1000000, 10000000};
    int nEnds[numIntervals] = {1000000, 10000000, 100000000};
    int nSteps[numIntervals] = {100000, 1000000, 10000000};
    gnuplot_ctrl *h1 = gnuplot_init();
    gnuplot_cmd(h1, "set terminal png size 2000,2000\n");
        char * setoutputcmd = malloc(sizeof *setoutputcmd * (strlen(outfile) + 20));
    sprintf(setoutputcmd, "set output \"%s\"\n", outfile);
    gnuplot_cmd(h1, setoutputcmd);
    gnuplot_cmd(h1, "set multiplot layout 2,2\n");
    for (int i = 0; i < numIntervals; ++i) {
        const int NPOINTS = (nEnds[i] - nStarts[i]) / nSteps[i] + 1;
        const int CPOINTS = (cEnd - cStart) / cStep + 1;

        bool(*built)[CPOINTS] = generate_FXLT_build_success_plot_helper(
            nStarts[i], nEnds[i], nSteps[i], cStart, cEnd, cStep, k);

        // print data to file
        FILE *temp = fopen("temp.dat", "w");
        // row header
        for (int j = 0; j < CPOINTS; ++j) {
            double c = cStart + j * cStep;
            fprintf(temp, ",c = %.3f", c);
        }
        // print matrix
        for (int x = 0; x < NPOINTS; ++x) {
            fprintf(temp, "\nn = %d", nStarts[i] + x * nSteps[i]);
            for (int y = 0; y < CPOINTS; ++y) {
                fprintf(temp, ", %d", built[x][y]);
            }
        }
        fflush(temp);
        fclose(temp);

        // plot on a heatmap
        gnuplot_cmd(h1, "set xlabel \"c\"\n");
        gnuplot_cmd(h1, "set ylabel \"n\"\n");
        gnuplot_cmd(h1,
                    "set title \"Build Time (s) vs. c (where m = c * n, n = "
                    "1,000,000)\"\n");

        gnuplot_cmd(h1, "set palette defined (0 'black', 1 'green')\n");
        // gnuplot_cmd(h1, "set cblabel \"build time (s)\"\n");
        // disable color bar
        gnuplot_cmd(h1, "unset colorbox\n");
        gnuplot_cmd(h1, "set nokey\n");
        gnuplot_cmd(h1, "set datafile separator comma\n");
        gnuplot_cmd(h1,
                    "plot \"temp.dat\" matrix rowheaders columnheaders using "
                    "1:2:3 with image\n");
    }
}

void testFXLT(int n, double c, int k, size_t flags) {
    SizedPointer *keys = malloc(sizeof *keys * n);
    void **vals = malloc(sizeof *vals * n);
    for (int i = 0; i < n; ++i) {
        char *intStr = malloc(sizeof *intStr * (i % 10 + 10));
        int len = sprintf(intStr, "%d long", i);
        keys[i].ptr = intStr;
        keys[i].size = len;
        vals[i] = (void *)(uintptr_t)((i & 1) + 1);
    }

    size_t kvmem = n * sizeof(void *);
    for (int j = 0; j < n; ++j) {
        kvmem += keys[j].size * sizeof(char);
    }

    struct timespec start, finish;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);
    struct fuseXORierLookupTable * fxlt = build_fuseXORierLT(n, keys, vals,
                                               c, k, FXLT_FLAG_PRINT_STATS | flags);

    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &finish);
    if (fxlt != NULL) {
        printf("FXLT with c=%f n=%d k=%d done in %fs\n", c, n, k,
               (finish.tv_sec - start.tv_sec) +
                   (double)(finish.tv_nsec - start.tv_nsec) / 1000000000.0);

    } else {
        printf("FXLT with c=%f n=%d k=%d failed\n", c, n, k);
        for (int i = 0; i < n; ++i) {
            free(keys[i].ptr);
        }
        free(keys);
        free(vals);
        return;
    }

    for (int i = 0; i < n; ++i) {
        uintptr_t v = fuseXORierLT_lookup(fxlt, keys[i]);
        if (v != (uintptr_t)vals[i]) {
            printf("Failed to lookup %s\n", keys[i].ptr);
            return;
        }
    }
    printf("All lookups successful\n");

    int fp = 0;
    for (int i = 0; i < n; ++i) {
        char *intStr = malloc(sizeof *intStr * (i % 10 + 10));
        int len = sprintf(intStr, "%d dif", i);
        SizedPointer sp = {.ptr = intStr, .size = len};
        uintptr_t v = fuseXORierLT_lookup(fxlt, sp);
        if (v != 0) {
            ++fp;
        }
    }
    printf("False Positive Rate: %f\n", (double)fp / n);

    delete_fuseXORierLT(fxlt);
    for (int i = 0; i < n; ++i) {
        free(keys[i].ptr);
    }
    free(keys);
    free(vals);
}

//------------------------------------------------------------------------------
// Main
//------------------------------------------------------------------------------

int main(void) {
    // generate_FXLT_vs_other_maps_plot("../../fig.png", 1000000, 10000000, 100000, 1.1, 4);
    // generate_spacial_coupling_advantage_plot("../../fig.png", 1.0, 1.4, .005, 10000000);
    // generate_FXLT_build_success_plot("../../fig.png", 1.02, 1.085, .005, 3);
    // generate_hash_advantage_plot("../../fig.png", 1000000, 6000000, 50000, 1.1, 4);
    
    printf("---Testing FXLT---\n");
    testFXLT(10000000, 1.075, 4, 0);
    printf("---Testing FXLT w/ Cache---\n");
    testFXLT(10000000, 1.075, 4, FXLT_FLAG_CACHE_HASHES);
    printf("---Testing FXLT k = 3---\n");
    testFXLT(10000000, 1.125, 3, 0);
    printf("---Testing FXLT k = 3 w/ Cache---\n");
    testFXLT(10000000, 1.125, 3, FXLT_FLAG_CACHE_HASHES);
    printf("---Testing FXLT w/ No Spatial Coupling---\n");
    testFXLT(10000000, 1.23, 3, FXLT_FLAG_NO_SPATIAL_COUPLING | FXLT_FLAG_ALLOW_RETRY_BUILD);

    return 0;
}
