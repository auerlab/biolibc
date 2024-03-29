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
#include "gff3.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seqid member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->seqid[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the seqid array
 *      new_seqid_element The new value for seqid[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char            new_seqid_element;
 *
 *      if ( bl_gff3_set_seqid_ae(&bl_gff, c, new_seqid_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_SEQID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_seqid_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char new_seqid_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->seqid[c] = new_seqid_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seqid member in a bl_gff3_t structure.
 *      Use this function to set seqid in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seqid to bl_gff3_ptr->seqid.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_seqid       The new value for seqid
 *      array_size      Size of the seqid array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char            new_seqid;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_seqid_cpy(&bl_gff, new_seqid, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_SEQID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_seqid_cpy(bl_gff3_t *bl_gff3_ptr, char new_seqid[], size_t array_size)

{
    if ( new_seqid == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->seqid, new_seqid, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of source member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->source[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the source array
 *      new_source_element The new value for source[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char            new_source_element;
 *
 *      if ( bl_gff3_set_source_ae(&bl_gff, c, new_source_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_SOURCE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_source_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char new_source_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->source[c] = new_source_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for source member in a bl_gff3_t structure.
 *      Use this function to set source in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_source to bl_gff3_ptr->source.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_source      The new value for source
 *      array_size      Size of the source array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char            new_source;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_source_cpy(&bl_gff, new_source, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_SOURCE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_source_cpy(bl_gff3_t *bl_gff3_ptr, char new_source[], size_t array_size)

{
    if ( new_source == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->source, new_source, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of type member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->type[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the type array
 *      new_type_element The new value for type[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char            new_type_element;
 *
 *      if ( bl_gff3_set_type_ae(&bl_gff, c, new_type_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_TYPE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_type_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char new_type_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->type[c] = new_type_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for type member in a bl_gff3_t structure.
 *      Use this function to set type in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_type to bl_gff3_ptr->type.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_type        The new value for type
 *      array_size      Size of the type array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char            new_type;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_type_cpy(&bl_gff, new_type, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_TYPE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_type_cpy(bl_gff3_t *bl_gff3_ptr, char new_type[], size_t array_size)

{
    if ( new_type == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->type, new_type, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for start member in a bl_gff3_t structure.
 *      Use this function to set start in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_start       The new value for start
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      int64_t         new_start;
 *
 *      if ( bl_gff3_set_start(&bl_gff, new_start)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_start(bl_gff3_t *bl_gff3_ptr, int64_t new_start)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->start = new_start;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for end member in a bl_gff3_t structure.
 *      Use this function to set end in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_end         The new value for end
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      int64_t         new_end;
 *
 *      if ( bl_gff3_set_end(&bl_gff, new_end)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_end(bl_gff3_t *bl_gff3_ptr, int64_t new_end)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->end = new_end;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for score member in a bl_gff3_t structure.
 *      Use this function to set score in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      score is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_score       The new value for score
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      double          new_score;
 *
 *      if ( bl_gff3_set_score(&bl_gff, new_score)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_score(bl_gff3_t *bl_gff3_ptr, double new_score)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->score = new_score;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for strand member in a bl_gff3_t structure.
 *      Use this function to set strand in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      strand is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_strand      The new value for strand
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char            new_strand;
 *
 *      if ( bl_gff3_set_strand(&bl_gff, new_strand)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_strand(bl_gff3_t *bl_gff3_ptr, char new_strand)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->strand = new_strand;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phase member in a bl_gff3_t structure.
 *      Use this function to set phase in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      phase is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_phase       The new value for phase
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char            new_phase;
 *
 *      if ( bl_gff3_set_phase(&bl_gff, new_phase)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_phase(bl_gff3_t *bl_gff3_ptr, char new_phase)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->phase = new_phase;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes member in a bl_gff3_t structure.
 *      Use this function to set attributes in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_attributes  The new value for attributes
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_attributes;
 *
 *      if ( bl_gff3_set_attributes(&bl_gff, new_attributes)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_attributes(bl_gff3_t *bl_gff3_ptr, char * new_attributes)

{
    if ( new_attributes == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->attributes = new_attributes;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of attributes member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->attributes[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the attributes array
 *      new_attributes_element The new value for attributes[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char *          new_attributes_element;
 *
 *      if ( bl_gff3_set_attributes_ae(&bl_gff, c, new_attributes_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_ATTRIBUTES_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_attributes_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char  new_attributes_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->attributes[c] = new_attributes_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes member in a bl_gff3_t structure.
 *      Use this function to set attributes in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_attributes to bl_gff3_ptr->attributes.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_attributes  The new value for attributes
 *      array_size      Size of the attributes array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_attributes;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_attributes_cpy(&bl_gff, new_attributes, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_ATTRIBUTES(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_attributes_cpy(bl_gff3_t *bl_gff3_ptr, char * new_attributes, size_t array_size)

{
    if ( new_attributes == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->attributes, new_attributes, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes_array_size member in a bl_gff3_t structure.
 *      Use this function to set attributes_array_size in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_attributes_array_size The new value for attributes_array_size
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          new_attributes_array_size;
 *
 *      if ( bl_gff3_set_attributes_array_size(&bl_gff, new_attributes_array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_attributes_array_size(bl_gff3_t *bl_gff3_ptr, size_t new_attributes_array_size)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->attributes_array_size = new_attributes_array_size;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes_len member in a bl_gff3_t structure.
 *      Use this function to set attributes_len in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_attributes_len The new value for attributes_len
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          new_attributes_len;
 *
 *      if ( bl_gff3_set_attributes_len(&bl_gff, new_attributes_len)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_attributes_len(bl_gff3_t *bl_gff3_ptr, size_t new_attributes_len)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->attributes_len = new_attributes_len;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_id member in a bl_gff3_t structure.
 *      Use this function to set feature_id in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_id is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_feature_id  The new value for feature_id
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_feature_id;
 *
 *      if ( bl_gff3_set_feature_id(&bl_gff, new_feature_id)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_id(bl_gff3_t *bl_gff3_ptr, char * new_feature_id)

{
    if ( new_feature_id == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->feature_id = new_feature_id;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_id member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->feature_id[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_id array
 *      new_feature_id_element The new value for feature_id[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_id_element;
 *
 *      if ( bl_gff3_set_feature_id_ae(&bl_gff, c, new_feature_id_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_FEATURE_ID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_id_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char  new_feature_id_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->feature_id[c] = new_feature_id_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_id member in a bl_gff3_t structure.
 *      Use this function to set feature_id in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_id to bl_gff3_ptr->feature_id.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_feature_id  The new value for feature_id
 *      array_size      Size of the feature_id array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_feature_id;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_feature_id_cpy(&bl_gff, new_feature_id, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_FEATURE_ID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_id_cpy(bl_gff3_t *bl_gff3_ptr, char * new_feature_id, size_t array_size)

{
    if ( new_feature_id == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->feature_id, new_feature_id, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_name member in a bl_gff3_t structure.
 *      Use this function to set feature_name in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_name is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_feature_name The new value for feature_name
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_feature_name;
 *
 *      if ( bl_gff3_set_feature_name(&bl_gff, new_feature_name)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_name(bl_gff3_t *bl_gff3_ptr, char * new_feature_name)

{
    if ( new_feature_name == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->feature_name = new_feature_name;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_name member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->feature_name[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_name array
 *      new_feature_name_element The new value for feature_name[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_name_element;
 *
 *      if ( bl_gff3_set_feature_name_ae(&bl_gff, c, new_feature_name_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_FEATURE_NAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_name_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char  new_feature_name_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->feature_name[c] = new_feature_name_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_name member in a bl_gff3_t structure.
 *      Use this function to set feature_name in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_name to bl_gff3_ptr->feature_name.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_feature_name The new value for feature_name
 *      array_size      Size of the feature_name array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_feature_name;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_feature_name_cpy(&bl_gff, new_feature_name, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_FEATURE_NAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_name_cpy(bl_gff3_t *bl_gff3_ptr, char * new_feature_name, size_t array_size)

{
    if ( new_feature_name == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->feature_name, new_feature_name, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_parent member in a bl_gff3_t structure.
 *      Use this function to set feature_parent in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_parent is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_feature_parent The new value for feature_parent
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_feature_parent;
 *
 *      if ( bl_gff3_set_feature_parent(&bl_gff, new_feature_parent)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_parent(bl_gff3_t *bl_gff3_ptr, char * new_feature_parent)

{
    if ( new_feature_parent == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->feature_parent = new_feature_parent;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_parent member in a bl_gff3_t
 *      structure. Use this function to set bl_gff3_ptr->feature_parent[c]
 *      in a bl_gff3_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_parent array
 *      new_feature_parent_element The new value for feature_parent[c]
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_parent_element;
 *
 *      if ( bl_gff3_set_feature_parent_ae(&bl_gff, c, new_feature_parent_element)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_FEATURE_PARENT_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_parent_ae(bl_gff3_t *bl_gff3_ptr, size_t c, char  new_feature_parent_element)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->feature_parent[c] = new_feature_parent_element;
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_parent member in a bl_gff3_t structure.
 *      Use this function to set feature_parent in a bl_gff3_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_parent to bl_gff3_ptr->feature_parent.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_feature_parent The new value for feature_parent
 *      array_size      Size of the feature_parent array.
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      char *          new_feature_parent;
 *      size_t          array_size;
 *
 *      if ( bl_gff3_set_feature_parent_cpy(&bl_gff, new_feature_parent, array_size)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF3_SET_FEATURE_PARENT(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_feature_parent_cpy(bl_gff3_t *bl_gff3_ptr, char * new_feature_parent, size_t array_size)

{
    if ( new_feature_parent == NULL )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff3_ptr->feature_parent, new_feature_parent, array_size);
	return BL_GFF3_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff3.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for file_pos member in a bl_gff3_t structure.
 *      Use this function to set file_pos in a bl_gff3_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      file_pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff3_ptr      Pointer to the structure to set
 *      new_file_pos    The new value for file_pos
 *
 *  Returns:
 *      BL_GFF3_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF3_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff3_t        bl_gff;
 *      long            new_file_pos;
 *
 *      if ( bl_gff3_set_file_pos(&bl_gff, new_file_pos)
 *              == BL_GFF3_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff3.h
 ***************************************************************************/

int     bl_gff3_set_file_pos(bl_gff3_t *bl_gff3_ptr, long new_file_pos)

{
    if ( false )
	return BL_GFF3_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff3_ptr->file_pos = new_file_pos;
	return BL_GFF3_DATA_OK;
    }
}
