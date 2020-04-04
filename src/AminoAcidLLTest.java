import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class AminoAcidLLTest {
 @Test

    public void testRNA() {
     String test1 = ("GUAAGUCGA*");
     AminoAcidLL list1 = AminoAcidLL.createFromRNASequence(test1);
    }

    @Test
    public void testAAL() {
     String test = ("GUAAGUCGA");
     AminoAcidLL test1 = AminoAcidLL.createFromRNASequence(test);
        char[] exp = {'V', 'S', 'R'};
     assertArrayEquals(exp, test1.aminoAcidList());
    }

    @Test
    public void testAAL2() {
        AminoAcidLL test2 = AminoAcidLL.createFromRNASequence("ACU");
        char[] exp = {'T'};
        assertArrayEquals(exp, test2.aminoAcidList());
    }
    @Test
    public void testCodonCompare() {
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCCAGGAGCAGG");
        AminoAcidLL.sort(test);
        int exp = 1;
        assertEquals(exp, test.codonCompare(test));
    }
    @Test
    public void testAminoAcidsort() {
     AminoAcidLL testSort = AminoAcidLL.createFromRNASequence("GGGAGGCUUGAG");
     char[] exp = {'E', 'G', 'L', 'R'};
     assertArrayEquals(exp, testSort.sort());
    }
    @Test
    public void testIsSorted() {
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCUACGGAGCUUCGGAGC");
        assertEquals(false,test.isSorted());
    }
    @Test
    public void testIsSorted1() {
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCGGACCACUCG");
        assertEquals(true,test.isSorted());
    }
    @Test
    public void testAminoAcidCompare() {
     AminoAcidLL list1 = AminoAcidLL.createFromRNASequence("GCCACGAGCAGG");
     AminoAcidLL.sort(list1);
     AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("GCCACGAGCAGC");
     AminoAcidLL.sort(list2);
     int exp = 1;
     assertEquals(exp, list1.aminoAcidCompare(list2));
    }
    @Test
    public void testAminoAcidCounts() {
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCCAGCAGCAGC");
        int [] exp = {1, 1, 2};
        assertArrayEquals(exp, test.aminoAcidCounts());
    }
    @Test
    public void testAminoAcidCounts2() {
        AminoAcidLL test = AminoAcidLL.createFromRNASequence("GCCACGAGCAGC");
        int [] exp = {1, 2, 1};
        assertArrayEquals(exp, test.aminoAcidCounts());
    }
}