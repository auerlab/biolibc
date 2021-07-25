/***************************************************************************
 *  This file is automatically generated by gen-get-set.  Be sure to keep
 *  track of any manual changes.
 *
 *  These generated functions are not expected to be perfect.  Check and
 *  edit as needed before adding to your code.
 ***************************************************************************/

#include <string.h>
#include <ctype.h>
#include "pos-list.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for array_size member in a bl_pos_list_t structure.
 *      Use this function to set array_size in a bl_pos_list_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      array_size is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_array_size is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the bl_bed_t structure to set
 *      new_array_size  The new value for array_size
 *
 *  Returns:
 *      BL_DATA_OK if the new value is acceptable and assigned
 *      BL_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      size_t          new_array_size;
 *
 *      if ( bl_pos_list_set_array_size(&bl_pos_list, new_array_size) == BL_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-25  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_array_size(bl_pos_list_t *bl_pos_list_ptr, size_t new_array_size)

{
    /* FIXME: Replace this with a proper sanity check */
    if ( 0 )
	return BL_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->array_size = new_array_size;
	return BL_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for count member in a bl_pos_list_t structure.
 *      Use this function to set count in a bl_pos_list_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      count is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_count is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the bl_bed_t structure to set
 *      new_count       The new value for count
 *
 *  Returns:
 *      BL_DATA_OK if the new value is acceptable and assigned
 *      BL_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      size_t          new_count;
 *
 *      if ( bl_pos_list_set_count(&bl_pos_list, new_count) == BL_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-25  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_count(bl_pos_list_t *bl_pos_list_ptr, size_t new_count)

{
    /* FIXME: Replace this with a proper sanity check */
    if ( 0 )
	return BL_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->count = new_count;
	return BL_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for positions member in a bl_pos_list_t structure.
 *      Use this function to set positions in a bl_pos_list_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      positions is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_positions is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the bl_bed_t structure to set
 *      new_positions   The new value for positions
 *
 *  Returns:
 *      BL_DATA_OK if the new value is acceptable and assigned
 *      BL_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      uint64_t *      new_positions;
 *
 *      if ( bl_pos_list_set_positions(&bl_pos_list, new_positions) == BL_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-25  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_positions(bl_pos_list_t *bl_pos_list_ptr, uint64_t * new_positions)

{
    /* FIXME: Replace this with a proper sanity check */
    if ( new_positions == NULL )
	return BL_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->positions = new_positions;
	return BL_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of positions member in a bl_pos_list_t
 *      structure. Use this function to set an element of the array
 *      positions in a bl_pos_list_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_POS_LIST_SET_POSITIONS_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_positions_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the bl_bed_t structure to set
 *      c               Subscript to the positions array
 *      new_positions_element The new value for positions[c]
 *
 *  Returns:
 *      BL_DATA_OK if the new value is acceptable and assigned
 *      BL_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      size_t          c;
 *      uint64_t *      new_positions_element;
 *
 *      if ( bl_pos_list_set_positions(&bl_pos_list, c, new_positions_element) == BL_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_POS_LIST_SET_POSITIONS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-25  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_positions_ae(bl_pos_list_t *bl_pos_list_ptr, size_t c, uint64_t  new_positions_element)

{
    /* FIXME: Replace this with a proper sanity check */
    if ( !isprint(new_positions_element) && (new_positions_element != '\0') )
	return BL_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->positions[c] = new_positions_element;
	return BL_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for positions member in a bl_pos_list_t structure.
 *      Use this function to set positions in a bl_pos_list_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_positions to ->positions.
 *
 *      Note that there is an equivalent macro BL_POS_LIST_SET_POSITIONS(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_positions is guaranteed by other means.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the bl_bed_t structure to set
 *      new_positions   The new value for positions
 *      array_size      Size of the positions array.
 *
 *  Returns:
 *      BL_DATA_OK if the new value is acceptable and assigned
 *      BL_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      uint64_t *      new_positions;
 *      size_t          array_size;
 *
 *      if ( bl_pos_list_set_positions(&bl_pos_list, new_positions, array_size) == BL_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_POS_LIST_SET_POSITIONS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-25  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_positions_cpy(bl_pos_list_t *bl_pos_list_ptr, uint64_t * new_positions, size_t array_size)

{
    /* FIXME: Replace this with a proper sanity check */
    if ( new_positions == NULL )
	return BL_DATA_OUT_OF_RANGE;
    else
    {
	{
	    size_t  c;
	    
	    // FIXME: Assuming all elements should be copied
	    for (c = 0; c < array_size; ++c)
		bl_pos_list_ptr->positions[c] = new_positions[c];
	}
	return BL_DATA_OK;
    }
}
