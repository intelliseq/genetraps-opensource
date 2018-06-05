package pl.intelliseq.explorare;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.junit.Test;

import pl.intelliseq.explorare.model.Results;
import pl.intelliseq.explorare.model.hpo.DiseaseGeneDictionary;
import pl.intelliseq.explorare.model.hpo.HpoOboParser;
import pl.intelliseq.explorare.model.hpo.HpoTerm;
import pl.intelliseq.explorare.model.hpo.HpoTree;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

public class DiseaseAutocomplete {

    @Test
    public void diseaseAutocomplete() throws IOException {

        DiseaseGeneDictionary diseaseGeneDictionary = new DiseaseGeneDictionary();
        HpoTree hpoTree = new HpoOboParser().getHpoTree();

        String firstLetters = "Ehlers";

        //Integer resultsCount = 30;


        Results results = new Results();

        for (String disease : diseaseGeneDictionary.getDiseases()) {
            if (disease.toLowerCase().startsWith(firstLetters.toLowerCase())) {
                results.addResult(disease);
                //if (results.sizeGreaterOrEqualTo(resultsCount))
                   // break;
            }
        }

        for (String disease : diseaseGeneDictionary.getDiseases()) {
            if (disease.toLowerCase().contains(firstLetters.toLowerCase())) {
                results.addResult(disease);
                //if (results.sizeGreaterOrEqualTo(resultsCount))
                 //   break;
            }
        }

        for (Object disease : results.getResults()) {
            System.out.println(disease);
        }


    }

}
