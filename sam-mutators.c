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
#include "sam.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of qname member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->qname[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the qname array
 *      new_qname_element The new value for qname[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char            new_qname_element;
 *
 *      if ( bl_sam_set_qname_ae(&bl_sam, c, new_qname_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QNAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qname_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_qname_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qname[c] = new_qname_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qname member in a bl_sam_t structure.
 *      Use this function to set qname in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_qname to bl_sam_ptr->qname.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qname       The new value for qname
 *      array_size      Size of the qname array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char            new_qname;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_qname_cpy(&bl_sam, new_qname, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QNAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qname_cpy(bl_sam_t *bl_sam_ptr, char new_qname[], size_t array_size)

{
    if ( new_qname == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->qname, new_qname, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for flag member in a bl_sam_t structure.
 *      Use this function to set flag in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      flag is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_flag        The new value for flag
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      unsigned        new_flag;
 *
 *      if ( bl_sam_set_flag(&bl_sam, new_flag)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_flag(bl_sam_t *bl_sam_ptr, unsigned new_flag)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->flag = new_flag;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of rname member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->rname[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the rname array
 *      new_rname_element The new value for rname[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char            new_rname_element;
 *
 *      if ( bl_sam_set_rname_ae(&bl_sam, c, new_rname_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rname_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_rname_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->rname[c] = new_rname_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for rname member in a bl_sam_t structure.
 *      Use this function to set rname in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_rname to bl_sam_ptr->rname.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_rname       The new value for rname
 *      array_size      Size of the rname array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char            new_rname;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_rname_cpy(&bl_sam, new_rname, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rname_cpy(bl_sam_t *bl_sam_ptr, char new_rname[], size_t array_size)

{
    if ( new_rname == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->rname, new_rname, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for pos member in a bl_sam_t structure.
 *      Use this function to set pos in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_pos         The new value for pos
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      int64_t         new_pos;
 *
 *      if ( bl_sam_set_pos(&bl_sam, new_pos)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_pos(bl_sam_t *bl_sam_ptr, int64_t new_pos)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->pos = new_pos;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for mapq member in a bl_sam_t structure.
 *      Use this function to set mapq in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      mapq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_mapq        The new value for mapq
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      unsigned char   new_mapq;
 *
 *      if ( bl_sam_set_mapq(&bl_sam, new_mapq)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_mapq(bl_sam_t *bl_sam_ptr, unsigned char new_mapq)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->mapq = new_mapq;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar member in a bl_sam_t structure.
 *      Use this function to set cigar in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      cigar is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar       The new value for cigar
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_cigar;
 *
 *      if ( bl_sam_set_cigar(&bl_sam, new_cigar)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar(bl_sam_t *bl_sam_ptr, char * new_cigar)

{
    if ( new_cigar == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar = new_cigar;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of cigar member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->cigar[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the cigar array
 *      new_cigar_element The new value for cigar[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char *          new_cigar_element;
 *
 *      if ( bl_sam_set_cigar_ae(&bl_sam, c, new_cigar_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_CIGAR_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_ae(bl_sam_t *bl_sam_ptr, size_t c, char  new_cigar_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar[c] = new_cigar_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar member in a bl_sam_t structure.
 *      Use this function to set cigar in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_cigar to bl_sam_ptr->cigar.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar       The new value for cigar
 *      array_size      Size of the cigar array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_cigar;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_cigar_cpy(&bl_sam, new_cigar, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_CIGAR(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_cpy(bl_sam_t *bl_sam_ptr, char * new_cigar, size_t array_size)

{
    if ( new_cigar == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->cigar, new_cigar, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of rnext member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->rnext[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the rnext array
 *      new_rnext_element The new value for rnext[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char            new_rnext_element;
 *
 *      if ( bl_sam_set_rnext_ae(&bl_sam, c, new_rnext_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNEXT_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rnext_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_rnext_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->rnext[c] = new_rnext_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for rnext member in a bl_sam_t structure.
 *      Use this function to set rnext in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_rnext to bl_sam_ptr->rnext.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_rnext       The new value for rnext
 *      array_size      Size of the rnext array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char            new_rnext;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_rnext_cpy(&bl_sam, new_rnext, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNEXT(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rnext_cpy(bl_sam_t *bl_sam_ptr, char new_rnext[], size_t array_size)

{
    if ( new_rnext == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->rnext, new_rnext, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for pnext member in a bl_sam_t structure.
 *      Use this function to set pnext in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      pnext is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_pnext       The new value for pnext
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      int64_t         new_pnext;
 *
 *      if ( bl_sam_set_pnext(&bl_sam, new_pnext)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_pnext(bl_sam_t *bl_sam_ptr, int64_t new_pnext)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->pnext = new_pnext;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for tlen member in a bl_sam_t structure.
 *      Use this function to set tlen in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      tlen is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_tlen        The new value for tlen
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      long            new_tlen;
 *
 *      if ( bl_sam_set_tlen(&bl_sam, new_tlen)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_tlen(bl_sam_t *bl_sam_ptr, long new_tlen)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->tlen = new_tlen;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_sam_t structure.
 *      Use this function to set seq in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq         The new value for seq
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_seq;
 *
 *      if ( bl_sam_set_seq(&bl_sam, new_seq)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq(bl_sam_t *bl_sam_ptr, char * new_seq)

{
    if ( new_seq == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq = new_seq;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seq member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->seq[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the seq array
 *      new_seq_element The new value for seq[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char *          new_seq_element;
 *
 *      if ( bl_sam_set_seq_ae(&bl_sam, c, new_seq_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_SEQ_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_ae(bl_sam_t *bl_sam_ptr, size_t c, char  new_seq_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq[c] = new_seq_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_sam_t structure.
 *      Use this function to set seq in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seq to bl_sam_ptr->seq.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq         The new value for seq
 *      array_size      Size of the seq array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_seq;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_seq_cpy(&bl_sam, new_seq, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_SEQ(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_cpy(bl_sam_t *bl_sam_ptr, char * new_seq, size_t array_size)

{
    if ( new_seq == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->seq, new_seq, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_sam_t structure.
 *      Use this function to set qual in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual        The new value for qual
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_qual;
 *
 *      if ( bl_sam_set_qual(&bl_sam, new_qual)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual(bl_sam_t *bl_sam_ptr, char * new_qual)

{
    if ( new_qual == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual = new_qual;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of qual member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->qual[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the qual array
 *      new_qual_element The new value for qual[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char *          new_qual_element;
 *
 *      if ( bl_sam_set_qual_ae(&bl_sam, c, new_qual_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QUAL_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_ae(bl_sam_t *bl_sam_ptr, size_t c, char  new_qual_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual[c] = new_qual_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_sam_t structure.
 *      Use this function to set qual in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_qual to bl_sam_ptr->qual.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual        The new value for qual
 *      array_size      Size of the qual array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_qual;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_qual_cpy(&bl_sam, new_qual, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QUAL(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_cpy(bl_sam_t *bl_sam_ptr, char * new_qual, size_t array_size)

{
    if ( new_qual == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->qual, new_qual, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar_array_size member in a bl_sam_t structure.
 *      Use this function to set cigar_array_size in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      cigar_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar_array_size The new value for cigar_array_size
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_cigar_array_size;
 *
 *      if ( bl_sam_set_cigar_array_size(&bl_sam, new_cigar_array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_array_size(bl_sam_t *bl_sam_ptr, size_t new_cigar_array_size)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar_array_size = new_cigar_array_size;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar_len member in a bl_sam_t structure.
 *      Use this function to set cigar_len in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      cigar_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar_len   The new value for cigar_len
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_cigar_len;
 *
 *      if ( bl_sam_set_cigar_len(&bl_sam, new_cigar_len)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_len(bl_sam_t *bl_sam_ptr, size_t new_cigar_len)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar_len = new_cigar_len;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_array_size member in a bl_sam_t structure.
 *      Use this function to set seq_array_size in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq_array_size The new value for seq_array_size
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_seq_array_size;
 *
 *      if ( bl_sam_set_seq_array_size(&bl_sam, new_seq_array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_array_size(bl_sam_t *bl_sam_ptr, size_t new_seq_array_size)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq_array_size = new_seq_array_size;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_len member in a bl_sam_t structure.
 *      Use this function to set seq_len in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq_len     The new value for seq_len
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_seq_len;
 *
 *      if ( bl_sam_set_seq_len(&bl_sam, new_seq_len)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_len(bl_sam_t *bl_sam_ptr, size_t new_seq_len)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq_len = new_seq_len;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual_array_size member in a bl_sam_t structure.
 *      Use this function to set qual_array_size in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual_array_size The new value for qual_array_size
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_qual_array_size;
 *
 *      if ( bl_sam_set_qual_array_size(&bl_sam, new_qual_array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_array_size(bl_sam_t *bl_sam_ptr, size_t new_qual_array_size)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual_array_size = new_qual_array_size;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual_len member in a bl_sam_t structure.
 *      Use this function to set qual_len in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual_len    The new value for qual_len
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_qual_len;
 *
 *      if ( bl_sam_set_qual_len(&bl_sam, new_qual_len)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_len(bl_sam_t *bl_sam_ptr, size_t new_qual_len)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual_len = new_qual_len;
	return BL_SAM_DATA_OK;
    }
}
