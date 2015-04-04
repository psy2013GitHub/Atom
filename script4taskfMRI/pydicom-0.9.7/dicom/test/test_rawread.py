# test_rawread.py
"""unittest tests for dicom.filereader module -- simple raw data elements"""
# Copyright (c) 2010 Darcy Mason
# This file is part of pydicom, relased under an MIT license.
#    See the file license.txt included with this distribution, also
#    available at http://pydicom.googlecode.com

from cStringIO import StringIO
import unittest
from dicom.filereader import data_element_generator
from dicom.values import convert_value
from dicom.sequence import Sequence

def hex2str(hexstr):
    """Return a bytestring rep of a string of hex rep of bytes separated by spaces"""
    return "".join((chr(int(x,16)) for x in hexstr.split()))

class RawReaderExplVRTests(unittest.TestCase):
    # See comments in data_element_generator -- summary of DICOM data element formats
    # Here we are trying to test all those variations

    def testExplVRLittleEndianLongLength(self):
        """Raw read: Explicit VR Little Endian long length......................"""
        # (0002,0001) OB 2-byte-reserved 4-byte-length, value 0x00 0x01
        infile = StringIO(hex2str("02 00 01 00 4f 42 00 00 02 00 00 00 00 01"))
        expected = ((2,1), 'OB', 2, '\00\01', 0xc, False, True)
        de_gen = data_element_generator(infile, is_implicit_VR=False, is_little_endian=True)
        got = de_gen.next()
        msg_loc = "in read of Explicit VR='OB' data element (long length format)"
        self.assertEqual(got, expected, "Expected: %r, got %r in %s" % (expected, got, msg_loc))
        # (0002,0002) OB 2-byte-reserved 4-byte-length, value 0x00 0x01
    def testExplVRLittleEndianShortLength(self):
        """Raw read: Explicit VR Little Endian short length....................."""
        # (0008,212a) IS 2-byte-length, value '1 '
        infile = StringIO(hex2str("08 00 2a 21 49 53 02 00 31 20"))
        expected = ((8,0x212a), 'IS', 2, '1 ', 0x8, False, True)
        de_gen = data_element_generator(infile, is_implicit_VR=False, is_little_endian=True)
        got = de_gen.next()
        msg_loc = "in read of Explicit VR='IS' data element (short length format)"
        self.assertEqual(got, expected, "Expected: %r, got %r in %s" % (expected, got, msg_loc))
    def testExplVRLittleEndianUndefLength(self):
        """Raw read: Expl VR Little Endian with undefined length................"""
        # (7fe0,0010), OB, 2-byte reserved, 4-byte-length (UNDEFINED)
        bytes1 = "e0 7f 10 00 4f 42 00 00 ff ff ff ff"
        bytes2 = " 41 42 43 44 45 46 47 48 49 4a"  # 'content'
        bytes3 = " fe ff dd e0 00 00 00 00"          # Sequence Delimiter
        bytes = bytes1 + bytes2 + bytes3
        infile = StringIO(hex2str(bytes))
        expected = ((0x7fe0,0x10), 'OB', 0xffffffffL, 'ABCDEFGHIJ', 0xc, False, True)
        de_gen = data_element_generator(infile, is_implicit_VR=False, is_little_endian=True)
        got = de_gen.next()
        msg_loc = "in read of undefined length Explicit VR ='OB' short value)"
        self.assertEqual(got, expected, "Expected: %r, got %r in %s" % (expected, got, msg_loc))

        # Test again such that delimiter crosses default 128-byte read "chunks", etc
        for multiplier in (116, 117, 118, 120):
            multiplier = 116
            bytes2b = bytes2 + " 00"*multiplier
            bytes = bytes1 + bytes2b + bytes3
            infile = StringIO(hex2str(bytes))
            expected = len('ABCDEFGHIJ'+'\0'*multiplier)
            de_gen = data_element_generator(infile, is_implicit_VR=False, is_little_endian=True)
            got = de_gen.next()
            got_len = len(got.value)
            msg_loc = "in read of undefined length Explicit VR ='OB' with 'multiplier' %d" % multiplier
            self.assertEqual(expected, got_len, "Expected value length %d, got %d in %s" % (expected, got_len, msg_loc))
            msg = "Unexpected value start with multiplier %d on Expl VR undefined length" % multiplier
            self.assert_(got.value.startswith('ABCDEFGHIJ\0'), msg)
            
class RawReaderImplVRTests(unittest.TestCase):
    # See comments in data_element_generator -- summary of DICOM data element formats
    # Here we are trying to test all those variations

    def testImplVRLittleEndian(self):
        """Raw read: Implicit VR Little Endian.................................."""
        # (0008,212a) {IS} 4-byte-length, value '1 '
        infile = StringIO(hex2str("08 00 2a 21 02 00 00 00 31 20"))
        expected = ((8,0x212a), None, 2, '1 ', 0x8, True, True)
        de_gen = data_element_generator(infile, is_implicit_VR=True, is_little_endian=True)
        got = de_gen.next()
        msg_loc = "in read of Implicit VR='IS' data element (short length format)"
        self.assertEqual(got, expected, "Expected: %r, got %r in %s" % (expected, got, msg_loc))
    def testImplVRLittleEndianUndefLength(self):
        """Raw read: Impl VR Little Endian with undefined length................"""
        # (7fe0,0010), OB, 2-byte reserved, 4-byte-length (UNDEFINED)
        bytes1 = "e0 7f 10 00 ff ff ff ff"
        bytes2 = " 41 42 43 44 45 46 47 48 49 4a"  # 'content'
        bytes3 = " fe ff dd e0 00 00 00 00"          # Sequence Delimiter
        bytes = bytes1 + bytes2 + bytes3
        infile = StringIO(hex2str(bytes))
        expected = ((0x7fe0,0x10), 'OW or OB', 0xffffffffL, 'ABCDEFGHIJ', 0x8, True, True)
        de_gen = data_element_generator(infile, is_implicit_VR=True, is_little_endian=True)
        got = de_gen.next()
        msg_loc = "in read of undefined length Implicit VR ='OB' short value)"
        self.assertEqual(got, expected, "Expected: %r, got %r in %s" % (expected, got, msg_loc))

        # Test again such that delimiter crosses default 128-byte read "chunks", etc
        for multiplier in (116, 117, 118, 120):
            multiplier = 116
            bytes2b = bytes2 + " 00"*multiplier
            bytes = bytes1 + bytes2b + bytes3
            infile = StringIO(hex2str(bytes))
            expected = len('ABCDEFGHIJ'+'\0'*multiplier)
            de_gen = data_element_generator(infile, is_implicit_VR=True, is_little_endian=True)
            got = de_gen.next()
            got_len = len(got.value)
            msg_loc = "in read of undefined length Implicit VR with 'multiplier' %d" % multiplier
            self.assertEqual(expected, got_len, "Expected value length %d, got %d in %s" % (expected, got_len, msg_loc))
            msg = "Unexpected value start with multiplier %d on Implicit VR undefined length" % multiplier
            self.assert_(got.value.startswith('ABCDEFGHIJ\0'), msg)

class RawSequenceTests(unittest.TestCase):
    # See DICOM standard PS3.5-2008 section 7.5 for sequence syntax
    def testEmptyItem(self):
        """Read sequence with a single empty item..............................."""
        # This is fix for issue 27
        bytes = (
             "08 00 32 10"    # (0008, 1032) SQ "Procedure Code Sequence"
            " 08 00 00 00"    # length 8
            " fe ff 00 e0"    # (fffe, e000) Item Tag
            " 00 00 00 00"    # length = 0
            ) + (             # --------------- end of Sequence
            " 08 00 3e 10"    # (0008, 103e) LO "Series Description"
            " 0c 00 00 00"    # length
            " 52 20 41 44 44 20 56 49 45 57 53 20"  # value
            )
        # "\x08\x00\x32\x10\x08\x00\x00\x00\xfe\xff\x00\xe0\x00\x00\x00\x00" # from issue 27, procedure code sequence (0008,1032)
        # bytes += "\x08\x00\x3e\x10\x0c\x00\x00\x00\x52\x20\x41\x44\x44\x20\x56\x49\x45\x57\x53\x20" # data element following
        
        fp = StringIO(hex2str(bytes))
        gen = data_element_generator(fp, is_implicit_VR=True, is_little_endian=True)
        raw_seq = gen.next()
        seq = convert_value("SQ", raw_seq)
        
        self.assert_(isinstance(seq, Sequence), "Did not get Sequence, got type '%s'" % type(seq))
        self.assert_(len(seq)==1, "Expected Sequence with single (empty) item, got %d item(s)" % len(seq))
        self.assert_(len(seq[0])==0, "Expected the sequence item (dataset) to be empty")
        elem2 = gen.next()
        self.assertEqual(elem2.tag, 0x0008103e, "Expected a data element after empty sequence item")
    
    def testImplVRLittleEndian_ExplicitLengthSeq(self):
        """Raw read: ImplVR Little Endian SQ with explicit lengths.............."""
        # Create a fictional sequence with bytes directly,
        #    similar to PS 3.5-2008 Table 7.5-1 p42
        bytes = (
            "0a 30 B0 00"    # (300a, 00b0) Beam Sequence
            " 40 00 00 00"    # length
                " fe ff 00 e0"    # (fffe, e000) Item Tag
                " 18 00 00 00"    # Item (dataset) Length
                " 0a 30 c0 00"    # (300A, 00C0) Beam Number
                " 02 00 00 00"    # length
                " 31 20"          # value '1 '
                " 0a 30 c2 00"    # (300A, 00C2) Beam Name
                " 06 00 00 00"    # length
                " 42 65 61 6d 20 31" # value 'Beam 1'
                # -------------
                " fe ff 00 e0"    # (fffe, e000) Item Tag
                " 18 00 00 00"    # Item (dataset) Length
                " 0a 30 c0 00"    # (300A, 00C0) Beam Number
                " 02 00 00 00"    # length
                " 32 20"          # value '2 '
                " 0a 30 c2 00"    # (300A, 00C2) Beam Name
                " 06 00 00 00"    # length
                " 42 65 61 6d 20 32" # value 'Beam 2'                
                )
                
        infile = StringIO(hex2str(bytes))
        de_gen = data_element_generator(infile, is_implicit_VR=True, is_little_endian=True)
        raw_seq = de_gen.next()
        seq = convert_value("SQ", raw_seq)

        # The sequence is parsed, but only into raw data elements. 
        # They will be converted when asked for. Check some:
        got = seq[0].BeamNumber
        self.assert_(got == 1, "Expected Beam Number 1, got %r" % got)
        got = seq[1].BeamName
        self.assert_(got == 'Beam 2', "Expected Beam Name 'Beam 2', got %s" % got)
        
    def testImplVRBigEndian_ExplicitLengthSeq(self):
        """Raw read: ImplVR BigEndian SQ with explicit lengths.................."""
        # Create a fictional sequence with bytes directly,
        #    similar to PS 3.5-2008 Table 7.5-1 p42
        bytes = (
            "30 0a 00 B0"    # (300a, 00b0) Beam Sequence
            " 00 00 00 40"    # length
                " ff fe e0 00"    # (fffe, e000) Item Tag
                " 00 00 00 18"    # Item (dataset) Length
                " 30 0a 00 c0"    # (300A, 00C0) Beam Number
                " 00 00 00 02"    # length
                " 31 20"          # value '1 '
                " 30 0a 00 c2"    # (300A, 00C2) Beam Name
                " 00 00 00 06"    # length
                " 42 65 61 6d 20 31" # value 'Beam 1'
                # -------------
                " ff fe e0 00"    # (fffe, e000) Item Tag
                " 00 00 00 18"    # Item (dataset) Length
                " 30 0a 00 c0"    # (300A, 00C0) Beam Number
                " 00 00 00 02"    # length
                " 32 20"          # value '2 '
                " 30 0a 00 c2"    # (300A, 00C2) Beam Name
                " 00 00 00 06"    # length
                " 42 65 61 6d 20 32" # value 'Beam 2'                
                )
                
        infile = StringIO(hex2str(bytes))
        de_gen = data_element_generator(infile, is_implicit_VR=True, is_little_endian=False)
        raw_seq = de_gen.next()
        seq = convert_value("SQ", raw_seq)

        # The sequence is parsed, but only into raw data elements. 
        # They will be converted when asked for. Check some:
        got = seq[0].BeamNumber
        self.assert_(got == 1, "Expected Beam Number 1, got %r" % got)
        got = seq[1].BeamName
        self.assert_(got == 'Beam 2', "Expected Beam Name 'Beam 2', got %s" % got)

    def testExplVRBigEndian_UndefinedLengthSeq(self):
        """Raw read: ExplVR BigEndian Undefined Length SQ......................."""
        # Create a fictional sequence with bytes directly,
        #    similar to PS 3.5-2008 Table 7.5-2 p42
        item1_value_bytes = "\1"*126
        item2_value_bytes = "\2"*222
        bytes = (
            "30 0a 00 B0"    # (300a, 00b0) Beam Sequence
            " 53 51"         # SQ
            " 00 00"         # reserved
            " ff ff ff ff"    # undefined length
                " ff fe e0 00"    # (fffe, e000) Item Tag
                " 00 00 00 18"    # Item (dataset) Length
                " 30 0a 00 c0"    # (300A, 00C0) Beam Number
                " 49 53"          # IS
                " 00 02"          # length
                " 31 20"          # value '1 '
                " 30 0a 00 c2"    # (300A, 00C2) Beam Name
                " 4c 4F"          # LO
                " 00 06"          # length
                " 42 65 61 6d 20 31" # value 'Beam 1'
                # -------------
                " ff fe e0 00"    # (fffe, e000) Item Tag
                " 00 00 00 18"    # Item (dataset) Length
                " 30 0a 00 c0"    # (300A, 00C0) Beam Number
                " 49 53"          # IS
                " 00 02"          # length
                " 32 20"          # value '2 '
                " 30 0a 00 c2"    # (300A, 00C2) Beam Name
                " 4C 4F"          # LO
                " 00 06"          # length
                " 42 65 61 6d 20 32" # value 'Beam 2'   
            " ff fe E0 dd"    # SQ delimiter
            " 00 00 00 00"    # zero length              
                )
                
        infile = StringIO(hex2str(bytes))
        de_gen = data_element_generator(infile, is_implicit_VR=False, is_little_endian=False)
        seq = de_gen.next()
        # Note seq itself is not a raw data element. 
        #     The parser does parse undefined length SQ

        # The sequence is parsed, but only into raw data elements. 
        # They will be converted when asked for. Check some:
        got = seq[0].BeamNumber
        self.assert_(got == 1, "Expected Beam Number 1, got %r" % got)
        got = seq[1].BeamName
        self.assert_(got == 'Beam 2', "Expected Beam Name 'Beam 2', got %s" % got)
        
if __name__ == "__main__":
    # import dicom
    # dicom.debug()
    unittest.main()
