import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

class AminoAcidLLTest {
 @Test

    public void testRNA() {
     String test1 = ("GUAAGUCGA*");
     AminoAcidLL list1 = AminoAcidLL.createFromRNASequence(test1);
    }

}