# filewriter.py
"""Write a dicom media file."""
# Copyright (c) 2008-2012 Darcy Mason
# This file is part of pydicom, released under a modified MIT license.
#    See the file license.txt included with this distribution, also
#    available at http://pydicom.googlecode.com

from struct import pack

import logging
logger = logging.getLogger('pydicom')

from dicom.UID import ExplicitVRLittleEndian, ImplicitVRLittleEndian, ExplicitVRBigEndian
from dicom.filebase import DicomFile
from dicom.datadict import dictionaryVR
from dicom.dataset import Dataset
from dicom.dataelem import DataElement
from dicom.tag import Tag, ItemTag, ItemDelimiterTag, SequenceDelimiterTag
from dicom.sequence import Sequence

def write_numbers(fp, data_element, struct_format):
    """Write a "value" of type struct_format from the dicom file.
    
    "Value" can be more than one number.
    
    struct_format -- the character format as used by the struct module.
    
    """
    endianChar = '><'[fp.is_little_endian]
    value = data_element.value
    if value == "":
        return  # don't need to write anything for empty string
    
    format_string = endianChar + struct_format
    try:
        try:
            value.append   # works only if list, not if string or number
        except: # is a single value - the usual case
            fp.write(pack(format_string, value))
        else:
            for val in value:
                fp.write(pack(format_string, val))
    except Exception, e:
        raise IOError, "%s\nfor data_element:\n%s" % (str(e), str(data_element))

def write_OBvalue(fp, data_element):
    """Write a data_element with VR of 'other byte' (OB)."""
    fp.write(data_element.value)

def write_OWvalue(fp, data_element):
    """Write a data_element with VR of 'other word' (OW).
    
    Note: This **does not currently do the byte swapping** for Endian state.
    
    """    
    # XXX for now just write the raw bytes without endian swapping
    fp.write(data_element.value)

def write_UI(fp, data_element):
    """Write a data_element with VR of 'unique identifier' (UI)."""    
    write_string(fp, data_element, '\0') # pad with 0-byte to even length

def multi_string(val):
    """Put a string together with delimiter if has more than one value"""
    if isinstance(val, (list, tuple)):
        return "\\".join(val)  # \ is escape chr, so "\\" gives single backslash
    else:
        return val

def write_string(fp, data_element, padding=' '):
    """Write a single or multivalued string."""
    val = multi_string(data_element.value)
    if len(val) % 2 != 0:
        val = val + padding   # pad to even length
    fp.write(val)

def write_number_string(fp, data_element, padding = ' '):
    """Handle IS or DS VR - write a number stored as a string of digits."""
    # If the DS or IS has an original_string attribute, use that, so that
    # unchanged data elements are written with exact string as when read from file 
    val = data_element.value   
    if isinstance(val, (list, tuple)):
        # XXX python >2.4 val = "\\".join([x.original_string if hasattr(x, 'original_string')
                        #  else str(x) for x in val])
        newval = []
        for x in val:
            if hasattr(x, 'original_string'):
                newval.append(x.original_string)
            else:
                newval.append(str(x))
        val = "\\".join(newval)
    else:
        # XXX python>2.4 val = val.original_string if hasattr(val, 'original_string') else str(val)
        if hasattr(val, 'original_string'):
            val = val.original_string
        else:
            val =  str(val)
    if len(val) % 2 != 0:
        val = val + padding   # pad to even length
    fp.write(val)
    
def write_data_element(fp, data_element):
    """Write the data_element to file fp according to dicom media storage rules."""
    fp.write_tag(data_element.tag)

    VR = data_element.VR
    if fp.is_explicit_VR:
        if len(VR) != 2:
            msg = "Cannot write ambiguous VR of '%s' for data element with tag %r." % (VR, data_element.tag)
            msg += "\nSet the correct VR before writing, or use an implicit VR transfer syntax"
            raise ValueError, msg
        fp.write(VR)
        if VR in ['OB', 'OW', 'OF', 'SQ', 'UT', 'UN']:
            fp.write_US(0)   # reserved 2 bytes
    if VR not in writers:
        raise NotImplementedError, "write_data_element: unknown Value Representation '%s'" % VR

    length_location = fp.tell() # save location for later.
    if fp.is_explicit_VR and VR not in ['OB', 'OW', 'OF', 'SQ', 'UT', 'UN']:
        fp.write_US(0)  # Explicit VR length field is only 2 bytes
    else:
        fp.write_UL(0xFFFFFFFFL)   # will fill in real length value later if not undefined length item
    
    try:
        writers[VR][0] # if writer is a tuple, then need to pass a number format
    except TypeError:
        writers[VR](fp, data_element) # call the function to write that kind of item
    else:
        writers[VR][0](fp, data_element, writers[VR][1])
    #  print DataElement(tag, VR, value)
    
    is_undefined_length = False
    if hasattr(data_element, "is_undefined_length") and data_element.is_undefined_length:
        is_undefined_length = True
    location = fp.tell()
    fp.seek(length_location)
    if fp.is_explicit_VR and VR not in ['OB', 'OW', 'OF', 'SQ', 'UT', 'UN']:
        fp.write_US(location - length_location - 2)  # 2 is length of US
    else:
        # write the proper length of the data_element back in the length slot, unless is SQ with undefined length.
        if not is_undefined_length:
            fp.write_UL(location - length_location - 4)  # 4 is length of UL
    fp.seek(location)  # ready for next data_element
    if is_undefined_length:
        fp.write_tag(SequenceDelimiterTag)
        fp.write_UL(0)  # 4-byte 'length' of delimiter data item
    
def write_dataset(fp, dataset):
    """Write a Dataset dictionary to the file. Return the total length written."""
    fpStart = fp.tell()
    # data_elements must be written in tag order
    tags = sorted(dataset.keys())
    for tag in tags:
        write_data_element(fp, dataset[tag])

    return fp.tell() - fpStart

def write_sequence(fp, data_element):
    """Write a dicom Sequence contained in data_element to the file fp."""
    # write_data_element has already written the VR='SQ' (if needed) and
    #    a placeholder for length"""
    sequence = data_element.value
    for dataset in sequence:
        write_sequence_item(fp, dataset)

def write_sequence_item(fp, dataset):
    """Write an item (dataset) in a dicom Sequence to the dicom file fp."""
    # see Dicom standard Part 5, p. 39 ('03 version)
    # This is similar to writing a data_element, but with a specific tag for Sequence Item
    fp.write_tag(ItemTag)   # marker for start of Sequence Item
    length_location = fp.tell() # save location for later.
    fp.write_UL(0xffffffffL)   # will fill in real value later if not undefined length
    write_dataset(fp, dataset)
    if getattr(dataset, "is_undefined_length_sequence_item", False):
        fp.write_tag(ItemDelimiterTag)
        fp.write_UL(0)  # 4-bytes 'length' field for delimiter item
    else: # we will be nice and set the lengths for the reader of this file
        location = fp.tell()
        fp.seek(length_location)
        fp.write_UL(location - length_location - 4)  # 4 is length of UL
        fp.seek(location)  # ready for next data_element

def write_UN(fp, data_element):
    """Write a byte string for an DataElement of value 'UN' (unknown)."""
    fp.write(data_element.value)

def write_ATvalue(fp, data_element):
    """Write a data_element tag to a file."""
    try:
        iter(data_element.value) # see if is multi-valued AT; # Note will fail if Tag ever derived from true tuple rather than being a long 
    except TypeError:
        tag = Tag(data_element.value)   # make sure is expressed as a Tag instance
        fp.write_tag(tag)
    else:
        tags = [Tag(tag) for tag in data_element.value]
        for tag in tags:
            fp.write_tag(tag)
            

def _write_file_meta_info(fp, meta_dataset):
    """Write the dicom group 2 dicom storage File Meta Information to the file.

    The file should already be positioned past the 128 byte preamble.
    Raises ValueError if the required data_elements (elements 2,3,0x10,0x12)
    are not in the dataset. If the dataset came from a file read with
    read_file(), then the required data_elements should already be there.
    """
    fp.write('DICM')

    # File meta info is always LittleEndian, Explicit VR. After will change these
    #    to the transfer syntax values set in the meta info
    fp.is_little_endian = True
    fp.is_explicit_VR = True

    if Tag((2,1)) not in meta_dataset:
        meta_dataset.add_new((2,1), 'OB', "\0\1")   # file meta information version
    
    # Now check that required meta info tags are present:
    missing = []
    for element in [2, 3, 0x10, 0x12]:
        if Tag((2, element)) not in meta_dataset:
            missing.append(Tag((2, element)))
    if missing:
        raise ValueError, "Missing required tags %s for file meta information" % str(missing)
    
    # Put in temp number for required group length, save current location to come back
    meta_dataset[(2,0)] = DataElement((2,0), 'UL', 0) # put 0 to start
    group_length_data_element_size = 12 # !based on DICOM std ExplVR
    group_length_tell = fp.tell()
    
    # Write the file meta datset, including temp group length
    length = write_dataset(fp, meta_dataset)
    group_length = length - group_length_data_element_size # counts from end of that
    
    # Save end of file meta to go back to
    end_of_file_meta = fp.tell()
    
    # Go back and write the actual group length
    fp.seek(group_length_tell)
    group_length_data_element = DataElement((2,0), 'UL', group_length)
    write_data_element(fp, group_length_data_element)
    
    # Return to end of file meta, ready to write remainder of the file
    fp.seek(end_of_file_meta)

def write_file(filename, dataset, WriteLikeOriginal=True):
    """Store a Dataset to the filename specified.
    
    Set dataset.preamble if you want something other than 128 0-bytes.
    If the dataset was read from an existing dicom file, then its preamble
    was stored at read time. It is up to you to ensure the preamble is still
    correct for its purposes.
    If there is no Transfer Syntax tag in the dataset,
       Set dataset.is_implicit_VR, and .is_little_endian
       to determine the transfer syntax used to write the file.
    WriteLikeOriginal -- True if want to preserve the following for each sequence 
        within this dataset:
        - preamble -- if no preamble in read file, than not used here
        - dataset.hasFileMeta -- if writer did not do file meta information,
            then don't write here either
        - seq.is_undefined_length -- if original had delimiters, write them now too,
            instead of the more sensible length characters
        - <dataset>.is_undefined_length_sequence_item -- for datasets that belong to a 
            sequence, write the undefined length delimiters if that is 
            what the original had
        Set WriteLikeOriginal = False to produce a "nicer" DICOM file for other readers,
            where all lengths are explicit.
    """

    # Decide whether to write DICOM preamble. Should always do so unless trying to mimic the original file read in
    preamble = getattr(dataset, "preamble", None) 
    if not preamble and not WriteLikeOriginal:
        preamble = "\0"*128
    
    file_meta = dataset.file_meta
    if file_meta is None:
        file_meta = Dataset()
    if 'TransferSyntaxUID' not in file_meta:
        if dataset.is_little_endian and dataset.is_implicit_VR:
            file_meta.add_new((2, 0x10), 'UI', ImplicitVRLittleEndian)
        elif dataset.is_little_endian and not dataset.is_implicit_VR:
            file_meta.add_new((2, 0x10), 'UI', ExplicitVRLittleEndian)
        elif dataset.is_big_endian and not dataset.is_implicit_VR:
            file_meta.add_new((2, 0x10), 'UI', ExplicitVRBigEndian)
        else:
            raise NotImplementedError, "pydicom has not been verified for Big Endian with Implicit VR"
        
    fp = DicomFile(filename,'wb')
    try:
        if preamble:
            fp.write(preamble)  # blank 128 byte preamble
            _write_file_meta_info(fp, file_meta) 
        
        # Set file VR, endian. MUST BE AFTER writing META INFO (which changes to Explict LittleEndian)
        fp.is_implicit_VR = dataset.is_implicit_VR
        fp.is_little_endian = dataset.is_little_endian
        
        write_dataset(fp, dataset)
    finally:
        fp.close()
        

WriteFile = write_file   # for backwards compatibility version <=0.9.2
writefile = write_file   # forgive user for missing underscore
        
# Map each VR to a function which can write it
# for write_numbers, the Writer maps to a tuple (function, struct_format)
#                                  (struct_format is python's struct module format)
writers = {'UL':(write_numbers,'L'), 'SL':(write_numbers,'l'),
           'US':(write_numbers,'H'), 'SS':(write_numbers, 'h'),
           'FL':(write_numbers,'f'), 'FD':(write_numbers, 'd'),
           'OF':(write_numbers,'f'),
           'OB':write_OBvalue, 'UI':write_UI,
           'SH':write_string,  'DA':write_string, 'TM': write_string,
           'CS':write_string,  'PN':write_string, 'LO': write_string,
           'IS':write_number_string,  'DS':write_number_string, 'AE': write_string,
           'AS':write_string,
           'LT':write_string,
           'SQ':write_sequence,
           'UN':write_UN,
           'AT':write_ATvalue,
           'ST':write_string,
           'OW':write_OWvalue,
           'US or SS':write_OWvalue,
           'OW/OB':write_OBvalue,
           'OB/OW':write_OBvalue,
           'OB or OW':write_OBvalue,
           'OW or OB':write_OBvalue,
           'DT':write_string,
           'UT':write_string,
           } # note OW/OB depends on other items, which we don't know at write time
