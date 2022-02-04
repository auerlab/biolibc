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
#include "gff.h"


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seqid member in a bl_gff_t
 *      structure. Use this function to set an element of the array
 *      seqid in a bl_gff_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_SEQID_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_seqid_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the seqid array
 *      new_seqid_element The new value for seqid[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char            new_seqid_element;
 *
 *      if ( bl_gff_set_seqid_ae(&bl_gff, c, new_seqid_element) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SEQID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_seqid_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_seqid_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->seqid[c] = new_seqid_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seqid member in a bl_gff_t structure.
 *      Use this function to set seqid in a bl_gff_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seqid to ->seqid.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_SEQID(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_seqid is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_seqid       The new value for seqid
 *      array_size      Size of the seqid array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_seqid;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_seqid_cpy(&bl_gff, new_seqid, array_size) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SEQID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_seqid_cpy(bl_gff_t *bl_gff_ptr, char new_seqid[], size_t array_size)

{
    if ( new_seqid == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->seqid, new_seqid, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of source member in a bl_gff_t
 *      structure. Use this function to set an element of the array
 *      source in a bl_gff_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_SOURCE_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_source_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the source array
 *      new_source_element The new value for source[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char            new_source_element;
 *
 *      if ( bl_gff_set_source_ae(&bl_gff, c, new_source_element) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SOURCE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_source_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_source_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->source[c] = new_source_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for source member in a bl_gff_t structure.
 *      Use this function to set source in a bl_gff_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_source to ->source.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_SOURCE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_source is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_source      The new value for source
 *      array_size      Size of the source array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_source;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_source_cpy(&bl_gff, new_source, array_size) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SOURCE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_source_cpy(bl_gff_t *bl_gff_ptr, char new_source[], size_t array_size)

{
    if ( new_source == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->source, new_source, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of type member in a bl_gff_t
 *      structure. Use this function to set an element of the array
 *      type in a bl_gff_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_TYPE_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_type_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the type array
 *      new_type_element The new value for type[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char            new_type_element;
 *
 *      if ( bl_gff_set_type_ae(&bl_gff, c, new_type_element) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_TYPE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_type_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_type_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->type[c] = new_type_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for type member in a bl_gff_t structure.
 *      Use this function to set type in a bl_gff_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_type to ->type.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_TYPE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_type is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_type        The new value for type
 *      array_size      Size of the type array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_type;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_type_cpy(&bl_gff, new_type, array_size) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_TYPE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_type_cpy(bl_gff_t *bl_gff_ptr, char new_type[], size_t array_size)

{
    if ( new_type == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->type, new_type, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for start member in a bl_gff_t structure.
 *      Use this function to set start in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      start is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_start is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_start       The new value for start
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      uint64_t        new_start;
 *
 *      if ( bl_gff_set_start(&bl_gff, new_start) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_start(bl_gff_t *bl_gff_ptr, uint64_t new_start)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->start = new_start;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for end member in a bl_gff_t structure.
 *      Use this function to set end in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      end is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_end is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_end         The new value for end
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      uint64_t        new_end;
 *
 *      if ( bl_gff_set_end(&bl_gff, new_end) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_end(bl_gff_t *bl_gff_ptr, uint64_t new_end)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->end = new_end;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for score member in a bl_gff_t structure.
 *      Use this function to set score in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      score is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_score is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_score       The new value for score
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      double          new_score;
 *
 *      if ( bl_gff_set_score(&bl_gff, new_score) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_score(bl_gff_t *bl_gff_ptr, double new_score)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->score = new_score;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for strand member in a bl_gff_t structure.
 *      Use this function to set strand in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      strand is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_strand is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_strand      The new value for strand
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_strand;
 *
 *      if ( bl_gff_set_strand(&bl_gff, new_strand) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_strand(bl_gff_t *bl_gff_ptr, char new_strand)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->strand = new_strand;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phase member in a bl_gff_t structure.
 *      Use this function to set phase in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      phase is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_phase is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_phase       The new value for phase
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_phase;
 *
 *      if ( bl_gff_set_phase(&bl_gff, new_phase) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_phase(bl_gff_t *bl_gff_ptr, char new_phase)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->phase = new_phase;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes member in a bl_gff_t structure.
 *      Use this function to set attributes in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_attributes is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_attributes  The new value for attributes
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_attributes;
 *
 *      if ( bl_gff_set_attributes(&bl_gff, new_attributes) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes(bl_gff_t *bl_gff_ptr, char * new_attributes)

{
    if ( new_attributes == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->attributes = new_attributes;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of attributes member in a bl_gff_t
 *      structure. Use this function to set an element of the array
 *      attributes in a bl_gff_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_ATTRIBUTES_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_attributes_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the attributes array
 *      new_attributes_element The new value for attributes[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_attributes_element;
 *
 *      if ( bl_gff_set_attributes_ae(&bl_gff, c, new_attributes_element) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_ATTRIBUTES_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_attributes_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->attributes[c] = new_attributes_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes member in a bl_gff_t structure.
 *      Use this function to set attributes in a bl_gff_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_attributes to ->attributes.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_ATTRIBUTES(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_attributes is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_attributes  The new value for attributes
 *      array_size      Size of the attributes array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_attributes;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_attributes_cpy(&bl_gff, new_attributes, array_size) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_ATTRIBUTES(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes_cpy(bl_gff_t *bl_gff_ptr, char * new_attributes, size_t array_size)

{
    if ( new_attributes == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->attributes, new_attributes, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_id member in a bl_gff_t structure.
 *      Use this function to set feature_id in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_id is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_feature_id is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_id  The new value for feature_id
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_id;
 *
 *      if ( bl_gff_set_feature_id(&bl_gff, new_feature_id) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_id(bl_gff_t *bl_gff_ptr, char * new_feature_id)

{
    if ( new_feature_id == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_id = new_feature_id;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_id member in a bl_gff_t
 *      structure. Use this function to set an element of the array
 *      feature_id in a bl_gff_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_FEATURE_ID_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_feature_id_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_id array
 *      new_feature_id_element The new value for feature_id[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_id_element;
 *
 *      if ( bl_gff_set_feature_id_ae(&bl_gff, c, new_feature_id_element) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_ID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_id_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_feature_id_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_id[c] = new_feature_id_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_id member in a bl_gff_t structure.
 *      Use this function to set feature_id in a bl_gff_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_id to ->feature_id.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_FEATURE_ID(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_feature_id is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_id  The new value for feature_id
 *      array_size      Size of the feature_id array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_id;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_feature_id_cpy(&bl_gff, new_feature_id, array_size) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_ID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_id_cpy(bl_gff_t *bl_gff_ptr, char * new_feature_id, size_t array_size)

{
    if ( new_feature_id == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->feature_id, new_feature_id, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for gene_name member in a bl_gff_t structure.
 *      Use this function to set gene_name in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      gene_name is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_gene_name is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_gene_name   The new value for gene_name
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_gene_name;
 *
 *      if ( bl_gff_set_gene_name(&bl_gff, new_gene_name) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_gene_name(bl_gff_t *bl_gff_ptr, char * new_gene_name)

{
    if ( new_gene_name == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->gene_name = new_gene_name;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of gene_name member in a bl_gff_t
 *      structure. Use this function to set an element of the array
 *      gene_name in a bl_gff_t variable from non-member functions.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_GENE_NAME_AE(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_gene_name_element is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the gene_name array
 *      new_gene_name_element The new value for gene_name[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_gene_name_element;
 *
 *      if ( bl_gff_set_gene_name_ae(&bl_gff, c, new_gene_name_element) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_GENE_NAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_gene_name_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_gene_name_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->gene_name[c] = new_gene_name_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for gene_name member in a bl_gff_t structure.
 *      Use this function to set gene_name in a bl_gff_t variable
 *      from non-member functions.  This function copies the array pointed to
 *      by new_gene_name to ->gene_name.
 *
 *      Note that there is an equivalent macro BL_GFF_SET_GENE_NAME(), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_gene_name is guaranteed by other means.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_gene_name   The new value for gene_name
 *      array_size      Size of the gene_name array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_gene_name;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_gene_name_cpy(&bl_gff, new_gene_name, array_size) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_GENE_NAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_gene_name_cpy(bl_gff_t *bl_gff_ptr, char * new_gene_name, size_t array_size)

{
    if ( new_gene_name == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->gene_name, new_gene_name, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for file_pos member in a bl_gff_t structure.
 *      Use this function to set file_pos in a bl_gff_t variable
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      file_pos is a pointer, data previously pointed to should
 *      generally be freed before calling this function to avoid memory
 *      leaks.
 *
 *      Note that there is an equivalent macro (), which performs
 *      this function with no data verification or function call overhead.
 *      Use the macro version to maximize performance where the validity
 *      of new_file_pos is guaranteed by other means.
 *      
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_file_pos    The new value for file_pos
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      long            new_file_pos;
 *
 *      if ( bl_gff_set_file_pos(&bl_gff, new_file_pos) == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-04  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_file_pos(bl_gff_t *bl_gff_ptr, long new_file_pos)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->file_pos = new_file_pos;
	return BL_GFF_DATA_OK;
    }
}
