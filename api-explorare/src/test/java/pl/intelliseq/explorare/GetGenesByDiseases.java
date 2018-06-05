package pl.intelliseq.explorare;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.Map.Entry;
import java.util.stream.Collectors;

import com.sun.scenario.effect.impl.sw.sse.SSEBlend_SRC_OUTPeer;
import org.junit.Test;

import org.omg.Messaging.SYNC_WITH_TRANSPORT;
import org.springframework.web.bind.annotation.RequestParam;
import pl.intelliseq.explorare.model.Results;
import pl.intelliseq.explorare.model.hpo.*;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

public class GetGenesByDiseases {

    @Test
    public void getGenesByDiseases() throws IOException {


        HpoOboParser hpoOboParser = new HpoOboParser();

        HpoTree hpoTree = hpoOboParser.getHpoTree();
        List<HpoTerm> hpoTerms = hpoOboParser.getTerms();

        List<String> finalGenesResult = new ArrayList<>();


        // Get Diseases list

        DiseaseGeneDictionary diseaseGeneDictionary = new DiseaseGeneDictionary();
        String firstLetters = "danlos";
        //String firstLetters = "ichthyosis";
        //String firstLetters = "netherton";
        //String firstLetters = "Larsson";
        //String firstLetters = "Arterial tortuosity syndrome";
        //String firstLetters = "blackfan";
        //String firstLetters = "marfan";



        //Integer resultsCount = 30;
        Results diseasesResults = new Results();

        for (String disease : diseaseGeneDictionary.getDiseases()) {
            if (disease.toLowerCase().startsWith(firstLetters.toLowerCase())) {
                diseasesResults.addResult(disease);
                //if (diseasesResults.sizeGreaterOrEqualTo(resultsCount))
                // break;
            }
        }

        for (String disease : diseaseGeneDictionary.getDiseases()) {
            if (disease.toLowerCase().contains(firstLetters.toLowerCase())) {
                diseasesResults.addResult(disease);
                //if (diseasesResults.sizeGreaterOrEqualTo(resultsCount))
                //   break;
            }
        }



        List <String> diseasesIds = new ArrayList<>();



        for (Object diseaseName : diseasesResults.getResults()) {

            // For each disease get phenotypes list

            //System.out.println("");
            //System.out.println("+ + + + + + + + + + + + + + + + + + + + + +");
            //System.out.println(diseaseName.toString());
            //System.out.println(diseaseGeneDictionary.getDiseaseByName(diseaseName.toString()).getId());
            //System.out.println(" ");


            System.out.println(diseaseGeneDictionary.getDiseaseByName(diseaseName.toString()).getId() + "\t" + diseaseName.toString());


            Results diseasePhenotypesResults = new Results();

            String diseaseId = diseaseGeneDictionary.getDiseaseByName(diseaseName.toString()).getId();


            //Set<HpoTerm> phenotypes = hpoTree.getPhenotypesByDisease(diseaseId.toString());
            Set<HpoTerm> phenotypes = makeSetSpecific(hpoTree.getPhenotypesByDisease(diseaseId.toString()));


            //for (HpoTerm term : phenotypes) {
            //    System.out.println(term.getId() + "\t" + term.getName());
            //}

            Map<String, Double> result = hpoTree.getDiseases(phenotypes)
                    .entrySet()
                    .stream()
                    .filter(map -> map.getValue() >= 0.3)
                    .collect(Collectors.toMap(p -> p.getKey(), p -> p.getValue()));

            // sorting
            result = result.entrySet().stream()
                    .sorted(Map.Entry.<String, Double>comparingByValue().reversed())
                    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue,
                            (e1, e2) -> e1, LinkedHashMap::new));


            Set<GeneResult> disiseasesGenesResult = new TreeSet<GeneResult>();


            for (Entry<String, Double> entry : result.entrySet()) {
                if (diseaseGeneDictionary.getDiseaseById(entry.getKey()) != null)

                    for (String gene : diseaseGeneDictionary.getDiseaseById(entry.getKey()).getGenes()) {
                        disiseasesGenesResult.add(new GeneResult(gene, 100 * entry.getValue()));
                        finalGenesResult.add(gene + "\t" + 100 * entry.getValue() + "\t" + diseaseName.toString());
                    }

            }

            // Print raw genes-list

            BufferedWriter outputWriter = null;
            outputWriter = new BufferedWriter(new FileWriter(firstLetters));

            for (String o : finalGenesResult) {
                outputWriter.write(o);
                outputWriter.newLine();
             }
            outputWriter.flush();
            outputWriter.close();


        }
    }

    private Set<HpoTerm> makeSetSpecific(Set<HpoTerm> hpoTerms) {

            Set<HpoTerm> out = new HashSet<HpoTerm>();

            for (HpoTerm term : hpoTerms) {

                boolean hasChildren = false;

                for (HpoTerm childTerm : hpoTerms) {
                    if (childTerm.isChildOf(term))
                        hasChildren = true;
                }

                if (!hasChildren)
                    out.add(term);
            }
            return out;

    }

}










