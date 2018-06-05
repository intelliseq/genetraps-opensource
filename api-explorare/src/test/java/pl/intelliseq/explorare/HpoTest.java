package pl.intelliseq.explorare;

import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.junit.Test;

import pl.intelliseq.explorare.model.hpo.HpoOboParser;
import pl.intelliseq.explorare.model.hpo.HpoTerm;
import pl.intelliseq.explorare.model.hpo.HpoTree;
import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

public class HpoTest {

	@Test
	public void parserTest() throws IOException {
		HpoTree hpoTree = new HpoOboParser().getHpoTree();
		assertTrue(hpoTree.getTerms().size() > 10000);
		
		assertTrue(hpoTree.getHpoTermById("HP:0003805").getGenes().contains("GNE"));
		
		
		Set <HpoTerm> osteogenesis = new HashSet<HpoTerm>();
		osteogenesis.add(hpoTree.getHpoTermById("HP:0002757"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0000592"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0012829"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0000365"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0000360"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0002321"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0012826"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0001382"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0002751"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0001084"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0002645"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0012825"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0000938"));
		osteogenesis.add(hpoTree.getHpoTermById("HP:0002980"));
		
		Map <String, Double> scores = hpoTree.getGenes(osteogenesis);
		
		int counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 30) break;
		}
		
		System.out.println();
		
		PhenoMarksParser parser = new PhenoMarksParser(hpoTree);
		String input = "";
		Set <HpoTerm> hpoTerms;
		
		System.out.println("Nonaka");
		
		input = "Argov and Yarom (1984) described a 'rimmed vacule myopathy' in Jews of Persian origin. The onset of this disorder usually occurred after the age of 20 years but before the middle of the fourth decade of life. Proximal and distal muscle weakness and wasting of the upper and lower limbs were progressive and resulted in severe incapacitation within 10 to 20 years. Despite this, there typically was sparing of the quadriceps muscles even in advanced stages of the disease, a feature unique to this form of inclusion body myopathy. It was not clear to the authors whether this disorder was primarily neurogenic or myopathic";
		hpoTerms = parser.tagInput(input).getHpoTerms();
		scores = hpoTree.getGenes(hpoTerms);

		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 15) break;
		}

		System.out.println();
		
		scores = hpoTree.getDiseases(hpoTerms);
		
		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 30) break;
		}
		
		System.out.println();
		
		System.out.println("Micro Warburg 1");
		
		input = "Warburg et al. (1993) used the designation Micro syndrome for an autosomal recessive syndrome comprising microcephaly, microcornea, congenital cataract, mental retardation, optic atrophy, and hypogenitalism. They described an affected brother and sister and their male cousin. The sibs were offspring of a consanguineous Pakistani marriage; the parents of the cousin denied consanguinity. Agenesis of the corpus callosum, prominent root of the nose, large anteverted ears, facial hypertrichosis, small pupils with posterior synechiae, hypotonia, mild to moderate spastic palsy with hip dislocations, and hormonal dysfunction, presumably of hypothalamic origin, were other features. The children were almost blind, whether or not the cataracts had been operated on. The electroretinographic responses indicated dysfunction of both retinal rods and cones, and the visual evoked potentials confirmed optic nerve atrophy. The children were late walkers and were incontinent of urine and stools. In the differential diagnosis, Warburg et al. (1993) considered COFS syndrome (214150), CAMAK/CAMFAK syndromes (212540), Martsolf syndrome (212720), lethal Rutledge syndrome (270400), and lethal Neu-Laxova syndrome (256520).  ";
		hpoTerms = parser.tagInput(input).getHpoTerms();
		scores = hpoTree.getGenes(hpoTerms);

		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 15) break;
		}

		System.out.println();
		
		scores = hpoTree.getDiseases(hpoTerms);
		
		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 30) break;
		}
		
		System.out.println();
		
		System.out.println("Osteogenesis");
		
		input = "Osteogenesis imperfecta (see Byers, 1993) is characterized chiefly by multiple bone fractures, usually resulting from minimal trauma. Affected individuals have blue sclerae, normal teeth, and normal or near-normal stature (for growth curves, see Vetter et al., 1992). Fractures are rare in the neonatal period; fracture tendency is constant from childhood to puberty, decreases thereafter, and often increases following menopause in women and after the sixth decade in men. Fractures heal rapidly with evidence of a good callus formation, and, with good orthopedic care, without deformity. Hearing loss of conductive or mixed type occurs in about 50% of families, beginning in the late teens and leading, gradually, to profound deafness, tinnitus, and vertigo by the end of the fourth to fifth decade. Additional clinical findings may be thin, easily bruised skin, moderate joint hypermobility and kyphoscoliosis, hernias, and arcus senilis. Mitral valve prolapse, aortic valvular insufficiency, and a slightly larger than normal aortic root diameter have been identified in some individuals (Hortop et al., 1986), but it is not clear that these disorders are significantly more frequent than in the general population. Radiologically, wormian bones are common but bone morphology is generally normal at birth, although mild osteopenia and femoral bowing may be present. Vertebral body morphology in the adult is normal initially, but often develops the classic 'cod-fish' appearance (Steinmann et al., 1991).";
		hpoTerms = parser.tagInput(input).getHpoTerms();
		scores = hpoTree.getGenes(hpoTerms);

		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 15) break;
		}

		System.out.println();
		
		scores = hpoTree.getDiseases(hpoTerms);
		
		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 30) break;
		}
		
		input = "Alzheimer disease is the most common form of progressive dementia in the elderly. It is a neurodegenerative disorder characterized by the neuropathologic findings of intracellular neurofibrillary tangles (NFT) and extracellular amyloid plaques that accumulate in vulnerable brain regions (Sennvik et al., 2000). Terry and Davies (1980) pointed out that the 'presenile' form, with onset before age 65, is identical to the most common form of late-onset or 'senile' dementia, and suggested the term 'senile dementia of the Alzheimer type' (SDAT).  Haines (1991) reviewed the genetics of AD. Selkoe (1996) reviewed the pathophysiology, chromosomal loci, and pathogenetic mechanisms of Alzheimer disease. Theuns and Van Broeckhoven (2000) reviewed the transcriptional regulation of the genes involved in Alzheimer disease.  ";
		hpoTerms = parser.tagInput(input).getHpoTerms();
		scores = hpoTree.getGenes(hpoTerms);

		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 15) break;
		}

		System.out.println();
		
		scores = hpoTree.getDiseases(hpoTerms);
		
		counter = 0;
		for (Entry <String, Double> score : scores.entrySet()) {
			System.out.println(score.getKey() + " " + score.getValue());
			if (counter++ > 30) break;
		}
		
	}

	
	
}
