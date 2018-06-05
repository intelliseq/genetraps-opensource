package pl.intelliseq.explorare;


import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.Map.Entry;
import org.junit.Test;

import pl.intelliseq.explorare.model.hpo.*;


public class GetGenesByPhenotypes {


    @Test
    public void parserTest() throws IOException {

        DecimalFormat df = new DecimalFormat("#.#");
        HpoTree hpoTree = new HpoOboParser().getHpoTree();
        DiseaseGeneDictionary diseaseGeneDictionary = new DiseaseGeneDictionary();


        // I. Add phenotypes //

        Set <HpoTerm> panel = new HashSet<HpoTerm>();

        //panel.add(hpoTree.getHpoTermById("HP:0001627"));
        //panel.add(hpoTree.getHpoTermById("HP:0001903"));
        //panel.add(hpoTree.getHpoTermById("HP:0000768"));
        //panel.add(hpoTree.getHpoTermById("HP:0000365"));
        panel.add(hpoTree.getHpoTermById("HP:0003549"));


        for (HpoTerm term: panel) {
            System.out.println(term.getId() + "\t" + term.getName());
        }



        Map <String, Double> scores = hpoTree.getGenes(panel);

        int count = 0;
        int count_30 = 0;

        for (Entry <String, Double> score : scores.entrySet()) {
            count ++;
            if (score.getValue() >= 0.3) {
                count_30 ++;
                System.out.println(score.getKey() + "\t" + df.format(score.getValue() * 100));

            }
        }

        System.out.println(count);
        System.out.println(count_30);



    }




}
