/***************************************************************************
 *  This file is automatically generated by gen-get-set.  Be sure to keep
 *  track of any manual changes.
 *
 *  These generated functions are not expected to be perfect.  Check and
 *  edit as needed before adding to your code.
 ***************************************************************************/

#include <string.h>
#include <ctype.h>
#include <stdbool.h>        // In case of bool
#include <stdint.h>         // In case of int64_t, etc
#include <xtend/string.h>   // strlcpy() on Linux
#include "overlap.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature1_len member in a bl_overlap_t structure.
 *      Use this function to set feature1_len in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature1_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_feature1_len The new value for feature1_len
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_feature1_len;
 *
 *      if ( bl_overlap_set_feature1_len(&bl_overlap, new_feature1_len)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_feature1_len(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_feature1_len
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->feature1_len = new_feature1_len;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature2_len member in a bl_overlap_t structure.
 *      Use this function to set feature2_len in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature2_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_feature2_len The new value for feature2_len
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_feature2_len;
 *
 *      if ( bl_overlap_set_feature2_len(&bl_overlap, new_feature2_len)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_feature2_len(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_feature2_len
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->feature2_len = new_feature2_len;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for overlap_start member in a bl_overlap_t structure.
 *      Use this function to set overlap_start in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      overlap_start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_overlap_start The new value for overlap_start
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_overlap_start;
 *
 *      if ( bl_overlap_set_overlap_start(&bl_overlap, new_overlap_start)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_overlap_start(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_overlap_start
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->overlap_start = new_overlap_start;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for overlap_end member in a bl_overlap_t structure.
 *      Use this function to set overlap_end in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      overlap_end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_overlap_end The new value for overlap_end
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_overlap_end;
 *
 *      if ( bl_overlap_set_overlap_end(&bl_overlap, new_overlap_end)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_overlap_end(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_overlap_end
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->overlap_end = new_overlap_end;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for overlap_len member in a bl_overlap_t structure.
 *      Use this function to set overlap_len in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      overlap_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_overlap_len The new value for overlap_len
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_overlap_len;
 *
 *      if ( bl_overlap_set_overlap_len(&bl_overlap, new_overlap_len)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_overlap_len(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_overlap_len
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->overlap_len = new_overlap_len;
	return BL_OVERLAP_DATA_OK;
    }
}
