class AminoAcidLL {
    char aminoAcid;
    String[] codons;
    int[] counts;
    AminoAcidLL next;

    AminoAcidLL() {

    }


    /********************************************************************************************/
    /* Creates a new node, with a given amino acid/codon
     * pair and increments the codon counter for that codon.
     * NOTE: Does not check for repeats!! */
    AminoAcidLL(String inCodon) {
        aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
        codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
        counts = new int[codons.length];
        incrCodons(inCodon);
        next = null;
    }


    private void incrCodons(String inCodon) {
        for (int i = 0; i < codons.length; i++) {
            if (codons[i].equals(inCodon)) {
                counts[i] += 1;
            }
        }
    }

    /********************************************************************************************/
    /* Recursive method that increments the count for a specific codon:
     * If it should be at this node, increments it and stops,
     * if not passes the task to the next node.
     * If there is no next node, add a new node to the list that would contain the codon.
     */
    private void addCodon(String inCodon) {
        if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
            incrCodons(inCodon);
        } else if (next != null) {
            next.addCodon(inCodon);
           // System.out.println("hello");
        } else {
            next = new AminoAcidLL(inCodon);
        }
    }

    /********************************************************************************************/
    /* Shortcut to find the total number of instances of this amino acid */
    private int totalCount() {
        int sum = 0;
        for (int i = 0; i < counts.length; i++) {
          sum += counts[i];
        }
        return sum;
    }

    /********************************************************************************************/
    /* helper method for finding the list difference on two matching nodes
     *  must be matching, but this is not tracked */
    private int totalDiff(AminoAcidLL inList) {
        return Math.abs(totalCount() - inList.totalCount());
    }


    /********************************************************************************************/
    /* helper method for finding the list difference on two matching nodes
     *  must be matching, but this is not tracked */
    private int codonDiff(AminoAcidLL inList) {
        int diff = 0;
        for (int i = 0; i < codons.length; i++) {
            diff += Math.abs(counts[i] - inList.counts[i]);
        }
        return diff;
    }

    /********************************************************************************************/
    /* Recursive method that finds the differences in **Amino Acid** counts.
     * the list *must* be sorted to use this method */
    public int aminoAcidCompare(AminoAcidLL inList) {
        //check to she if it is sorted
        if(!inList.isSorted()) {
            inList.sort(inList);
        }
        //Checks to see if it is at the end of both lists
        if (inList.next == null && next == null) {
            return codonDiff(inList);
        }
        //Check to see if it is at the end of the Amino Acid
        if (next == null) {
            return inList.totalCount() + aminoAcidCompare(inList);
        }

        if (inList.next == null) {
            return totalCount() + next.aminoAcidCompare(inList);
        }
        return totalDiff(inList) + next.aminoAcidCompare(inList.next);
    }

    /********************************************************************************************/
    /* Same ad above, but counts the codon usage differences
     * Must be sorted. */
    public int codonCompare(AminoAcidLL inList) {
        //Sends to sort if not sorted
        if(!inList.isSorted()){
            inList.sort(inList);
        }
        //Base case after going through the codons
        if(inList.next == null && next == null) {
            return codonDiff(inList);
        }
        //goes though the starting list
        if(next == null) {
            return inList.totalCount() + codonCompare(inList);
        }
        //goes through the list
        if(inList.next == null) {
            return totalCount() + codonCompare(inList);
        }
        return totalDiff(inList) + next.codonCompare(inList);
    }

    /********************************************************************************************/
    /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
    public char[] aminoAcidList() {
        char[] ret;
        //Checks to see if the list is complete
        if (next == null) {
            return new char[]{aminoAcid};
        } else {
            //Adds to a new array if the list is not completed
            char[] a = next.aminoAcidList();
            ret = new char[a.length + 1];
            ret[0] = aminoAcid;
            for (int i = 1; i < ret.length - 1; i++) {
                ret[i] = a[i - 1];
            }
        }
        return ret;
    }

    /********************************************************************************************/
    /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
    public int[] aminoAcidCounts() {
        int[] total;
        //this is the base case
        if(next == null) {
            return new int[totalCount()];
        }
        else {
            //Adds to the counter
            int[] a = next.aminoAcidCounts();
            total = new int[a.length + 1];
            total[0] = totalCount();
            for (int i = 1; i < total.length; i++) {
                total[i] = a[i];
            }
        }
        return total;
    }


    /********************************************************************************************/
    /* recursively determines if a linked list is sorted or not */
    public boolean isSorted() {
        if(next == null) {
            return true;
        }
        if(aminoAcid > next.aminoAcid) {
            return false;
        }

        return next.isSorted();
    }


    /********************************************************************************************/
    /* Static method for generating a linked list from an RNA sequence */
    public static AminoAcidLL createFromRNASequence(String inSequence) {
        AminoAcidLL head = new AminoAcidLL();
        boolean test = true;
        //Checks if there is a star and breaks the loop
        if (inSequence.substring(0, 3).charAt(0) == '*') {
            head.addCodon(inSequence.substring(0, 3));
            test = false;
        } else {
            head.addCodon(inSequence.substring(0, 3));
        }
            for (int i = 3; i < inSequence.length() - 2 && test; i += 3) {
                System.out.println(inSequence.substring(i, i + 3));
                if (inSequence.charAt(i) == '*') {
                    test = false;
                } else {
                    head.addCodon(inSequence.substring(i, i + 3));
                }
            }
        return head;
    }
//str.substring(0, 2);

    /********************************************************************************************/
    /* sorts a list by amino acid character*/
    public static AminoAcidLL sort(AminoAcidLL inList) {
        //Checks if the list is empty
        if(inList == null || inList.next == null) {
            return inList;
        }
        char temp;
        //selection sort method
        for(AminoAcidLL i = inList; i.next != null; i = i.next) {
            for(AminoAcidLL j = i.next; j != null; j = j.next) {
                temp = i.aminoAcid;
                i.aminoAcid = j.aminoAcid;
                j.aminoAcid = temp;
            }
        }

        return inList;
    }
    }
