#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/**
 * @brief enum of the color of the node
 */
enum color {
    red,
    black
};

enum insertFixupCase {
    root,
    uncleRed,
    triangle,
    line,
    good
};
enum insertFixupOrientation {
    case_left,
    case_right
};
struct insertFixupInfo {
    enum insertFixupCase fc;
    enum insertFixupOrientation fo;
};

/**
 * @brief a class that represents a rbtree node
 */
struct rbnode {
    /**
     * @brief a pointer to the pointer of the tree object's root
     */
    struct rbnode **root; // = nullptr;
    /**
     * @brief this nodes color
     */
    enum color c; // = red;
    /**
     * @brief the left child of this node
     */
    struct rbnode *l; // = nullptr;
    /**
     * @brief the right child of this node
     */
    struct rbnode *r; // = nullptr;
    /**
     * @brief the parent of this node
     */
    struct rbnode *p; // = nullptr;
    /**
     * @brief the key of this node
     */
    void * key;
    /**
     * @brief the value stored in this node
     */
    void * value;

    struct rbnode_VTABLE_s *vtable;
};

struct rbnode_VTABLE_s {
    /**
     * @brief provides a total ordering of keys.
     * K must overload operator-, ie. of a numeric type or IComparable
     *
     * @param[in] other an object of type K to compare
     * @return int 0 if equal, positive if this is greater, negative if other is greater
     */
    int (*compare)(struct rbnode *this, void * other);
    /**
     * @brief delete this node from the tree, fixing pointers accordingly
     *
     * @param[in] parentsChild pointer to the pointer the parent holds of the node to delete (used to update the parent's child pointer easily)
     */
    void (*commitApoptosis)(struct rbnode *this, struct rbnode **parentsChild);
    /**
     * @brief Get this node's uncle, if it exists
     *
     * @return struct rbnode* the uncle
     */
    struct rbnode *(*getUncle)(struct rbnode *this);
    /**
     * @brief Get the insert case
     *
     * @return insertFixupInfo the fixup case info
     */
    struct insertFixupInfo (*getInsertCase)(struct rbnode *this);
    /**
     * @brief the black magic to fix the red black tree after an insertion
     *
     */
    void (*insertFixup)(struct rbnode *this);
    /**
     * @brief perform a left rotation around this node
     */
    void (*rotateLeft)(struct rbnode *this);
    /**
     * @brief perform a right rotation around this node
     */
    void (*rotateRight)(struct rbnode *this);
    /**
     * @brief checks if a node is black (nil nodes are black)
     *
     * @param[in] n the node to check
     * @return true if it is black
     * @return false otherwise
     */
    int (*static_isBlack)(struct rbnode *n);
    /**
     * @brief checks if a node is red (nil nodes are black)
     *
     * @param[in] n the node to check
     * @return true if it is red
     * @return false otherwise
     */
    int (*static_isRed)(struct rbnode *n);
    /**
     * @brief Get the dbnode's sibling, if it exists
     *
     * @param[in] dbnode the double black node that needs fixing
     * @param[in] parent the parent of the dbnode (needed because the db node could be a nil node)
     * @return struct rbnode* the sibling
     */
    struct rbnode *(*static_getSibling)(struct rbnode *dbnode, struct rbnode *parent);
    /**
     * @brief Get the dbnode's sibling's child furthest from the dbnode, if it exists
     *
     * @param[in] dbnode the double black node that needs fixing
     * @param[in] parent the parent of the dbnode (needed because the db node could be a nil node)
     * @return struct rbnode* the sibling's far child
     */
    struct rbnode *(*static_getSiblingFarChild)(struct rbnode *dbnode, struct rbnode *parent);
    /**
     * @brief Get the dbnode's sibling's child nearest to the dbnode, if it exists
     *
     * @param[in] dbnode the double black node that needs fixing
     * @param[in] parent the parent of the dbnode (needed because the db node could be a nil node)
     * @return struct rbnode* the sibling's near child
     */
    struct rbnode *(*static_getSiblingNearChild)(struct rbnode *dbnode, struct rbnode *parent);
    /**
     * @brief the black magic to fix the red black tree after a removal
     * cases and how to fix them referenced from here: https://medium.com/analytics-vidhya/deletion-in-red-black-rb-tree-92301e1474ea
     *
     * @param[in] dbnode the double black node that needs fixing
     * @param[in] parent the parent of the dbnode (needed because the db node could be a nil node)
     */
    void (*static_doubleBlackFixup)(struct rbnode *dbnode, struct rbnode *parent);
    /**
     * @brief update the root reference of the tree object
     * @param[in] root a pointer to the new root
     */
    void (*setRoot)(struct rbnode *this, struct rbnode *root);

    //    public:
    /**
     * @brief Construct a new rbnode object
     *
     * @param[in] key it's key
     * @param[in] value it's value
     * @param[in] parent it's parent
     * @param[in] r a pointer to the root node pointer, in case it needs to be updated
     */
    /**
     * @brief print the tree rooted at this node using a reverse pre order traversal
     *
     * @param[in] depth the depth of this node
     */
    void (*printTree)(struct rbnode *this, int depth, FILE * out, void (*printKey)(void * thing,FILE * out),void (*printVal)(void * thing,FILE * out));

    int (*getKeys)(struct rbnode *this, void ***keys);
    /**
     * @brief makes a copy of the tree rooted at this node
     *
     * @param[in] parent the pointer to the parent of this copy
     * @return rbnode* the copy of the tree rooted here
     */
    struct rbnode *(*copy)(struct rbnode *this, struct rbnode *parent);
    /**
     * @brief check if the node with this key is in the tree rooted at this node
     *
     * @param[in] key the key to look for
     * @return true if the key is in this tree
     * @return false otherwise
     */
    int (*contains)(struct rbnode *this, void * key);
    /**
     * @brief get the value of the node with the key in the tree rooted at this node
     *
     * NOTE: this method has undefined behavior if the key is not in this subtree
     *
     * @param[in] key the key to look for
     * @return V the value of the key if it exists in the tree, undefined value otherwise
     */
    void * (*getValue)(struct rbnode *this, void * key);
    /**
     * @brief gets the key's value if it is in the tree rooted at this node
     *
     * @param[in] key the key to look for
     * @param[out] out the location to put the value if found
     * @return true if the key is in this tree
     * @return false otherwise
     */
    int (*tryGetValue)(struct rbnode *this, void * key, void * *out);
    /**
     * @brief set the value of a node with corrsponding key in the tree rooted at this node
     *
     * @param[in] key the key of the node to set
     * @param[in] value the value to set
     * @return true if the node was in the tree and set
     * @return false otherwise
     */
    int (*setValue)(struct rbnode *this, void * key, void * value);
    /**
     * @brief insert a key value pair into the tree rooted at this node
     *
     * @param[in] key the key of the node to insert
     * @param[in] value the value to be stored in the node
     * @return true if the node wasn't in the subtree and thus inserted
     * @return false otherwise
     */
    int (*insert)(struct rbnode *this, void * key, void * value);
    /**
     * @brief remove a node with a key from the tree rooted at this node
     *
     * @param[in] key the key of the node to remove
     * @param[in] parentsChild a pointer to the pointer the parent holds of the node to delete (used to update the parent's child pointer easily)
     * @return true if the node was in the subtree and deleted
     * @return false otherwise
     */
    int (*remove)(struct rbnode *this, void * key, struct rbnode **parentsChild);
    /**
     * @brief get the black height of the tree rooted at this node
     *
     * @return int the bh if it is valid or -1 if the black height property is broken
     */
    int (*bh)(struct rbnode *this);
    /**
     * @brief check to make sure the red property holds in the tree rooted at this node
     * true in the tree rooted at this node
     *
     * @return true if it does
     * @return false otherwise
     */
    int (*rp)(struct rbnode *this);
};

void destroy_rbnode(struct rbnode *this);
struct rbnode *construct_rbnode(void * key, void * value, struct rbnode *parent, struct rbnode **r);