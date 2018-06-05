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

public class GetGenesByDiseasesIds {

    @Test
    public void getGenesByDiseasesIds() throws IOException {

        HpoOboParser hpoOboParser = new HpoOboParser();

        HpoTree hpoTree = hpoOboParser.getHpoTree();
        List<HpoTerm> hpoTerms = hpoOboParser.getTerms();

        List<String> finalGenesResult = new ArrayList<>();



        DiseaseGeneDictionary diseaseGeneDictionary = new DiseaseGeneDictionary();


        // Get Diseases list
        Results diseasesResults = new Results();

        List<String> diseasesIds = new ArrayList<>();
        diseasesIds.add("OMIM:113610");
        diseasesIds.add("OMIM:159580");
        diseasesIds.add("OMIM:182610");
        diseasesIds.add("OMIM:182950");
        diseasesIds.add("OMIM:231000");
        diseasesIds.add("OMIM:238970");
        diseasesIds.add("OMIM:248900");
        diseasesIds.add("OMIM:260600");
        diseasesIds.add("OMIM:275900");
        diseasesIds.add("OMIM:277580");
        diseasesIds.add("OMIM:300166");
        diseasesIds.add("OMIM:300894");
        diseasesIds.add("OMIM:312910");
        diseasesIds.add("OMIM:312920");
        diseasesIds.add("OMIM:608804");
        diseasesIds.add("OMIM:609136");
        diseasesIds.add("OMIM:612319");
        diseasesIds.add("OMIM:613280");
        diseasesIds.add("OMIM:613672");
        diseasesIds.add("OMIM:614487");
        diseasesIds.add("OMIM:614877");
        diseasesIds.add("OMIM:615157");
        diseasesIds.add("OMIM:615643");
        diseasesIds.add("OMIM:616053");
        diseasesIds.add("OMIM:616680");
        diseasesIds.add("ORPHA:2710");
        diseasesIds.add("ORPHA:2815");
        diseasesIds.add("ORPHA:397725");
        diseasesIds.add("ORPHA:67047");
        diseasesIds.add("ORPHA:726");
        diseasesIds.add("ORPHA:93473");
        diseasesIds.add("ORPHA:93474");

         for (String diseaseIds : diseasesIds ){

             System.out.println(diseaseGeneDictionary.getDiseaseById(diseaseIds).getPrimaryName());
             for (String gene : diseaseGeneDictionary.getDiseaseById(diseaseIds).getGenes()) {
                 System.out.println(gene + "\t");
             }

         }


      //  Set<HpoTerm> phenotypes = new HashSet<>();
      //  phenotypes.add(hpoTree.getHpoTermById("HP:0002313"));


      //  Map<String, Double> result = hpoTree.getGenes(phenotypes)
       //         .entrySet()
       //         .stream()
        //        .filter(map -> map.getValue() >= 0.3)
       ///         .collect(Collectors.toMap(p -> p.getKey(), p -> p.getValue()));


       /// for (Entry<String, Double> entry : result.entrySet()) {
///
       //             System.out.println(entry.getKey()+ "\t" + 100 * entry.getValue());

       // }




/*

        for (String diseaseId : diseasesIds) {

            // For each disease get phenotypes list

            System.out.println("");
            System.out.println("+ + + + + + + + + + + + + + + + + + + + + +");
            System.out.println(diseaseGeneDictionary.getDiseaseById(diseaseId).getPrimaryName());
            System.out.println(diseaseId);
            System.out.println(" ");


            //System.out.println(diseaseGeneDictionary.getDiseaseByName(diseaseName.toString()).getId() + "\t" + diseaseName.toString());


            Results diseasePhenotypesResults = new Results();


            //Set<HpoTerm> phenotypes = hpoTree.getPhenotypesByDisease(diseaseId.toString());
            //Set<HpoTerm> phenotypes = makeSetSpecific(hpoTree.getPhenotypesByDisease(diseaseId));


            //for (HpoTerm term : phenotypes) {
            //    System.out.println(term.getId() + "\t" + term.getName());
            //}

/*
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
            */
            // Print raw genes-list

            //BufferedWriter outputWriter = null;
           // outputWriter = new BufferedWriter(new FileWriter(firstLetters));

            //for (String o : finalGenesResult) {
            //    outputWriter.write(o);
            //    outputWriter.newLine();
            // }
           // outputWriter.flush();
           // outputWriter.close();



        }

        private Set<HpoTerm> makeSetSpecific (Set < HpoTerm > hpoTerms) {

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










