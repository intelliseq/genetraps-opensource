package pl.intelliseq.explorare;

import static org.junit.Assert.assertTrue;

import java.io.IOException;

import org.junit.Test;

import pl.intelliseq.explorare.model.phenoMarks.PhenoMarksParser;

public class HpoOboParserTest {

	@Test
	public void parserTest() throws IOException {
		String input = ""; 
		PhenoMarksParser parser = new PhenoMarksParser();
		
		input = "The patient is 15 y.o boy with no prior somatic problems. The first prodromal symptoms occurred in June, 2015. The patient\n\r had transient no provoked pain in lower extremities and in right lower quadrant of the abdomen. It resolved spontaneously. Then in August, 2015 a severe pain occurred in the back and again in the right lower quadrant of the abdomen. This time the patient was transferred to ER. ECHO showed dissection of the aorta type DeBackeyIII/Stanford B. The valves and Ascending aorta appeared normal. CT of the abdomen showed the dissection started 3cm beneath of the origin of the Left subclavian artery and ended 27 cm into the right common iliac artery; infarct of the lower pole of the right kidney; double supply of the right kidney –lower pole from the false channel of the aorta and the upper pole from the true channel of the aorta; left iliac common artery is supplied from the true channel. As so far, the patient is treated symptomatically. In addition, he has been treated for the past 10 years for obsessive-compulsive and ADHD type of behavior. Medications: perazine, fluoxetine, haroperidol (on demand). Physical examination beside macrocephaly, mild scoliosis, flat foot, and abdominal tenderness, was normal.";
		assertTrue(parser.tagInput(input).getMarkedText().contains("HP:0002647"));
		//System.out.println(parser.tagInput(input).getMarkedText());
		
		//input = "Warburg et al. (1993) used the designation Micro syndrome for an autosomal recessive syndrome comprising microcephaly, microcornea, congenital cataract, mental retardation, optic atrophy, and hypogenitalism. They described an affected brother and sister and their male cousin. The sibs were offspring of a consanguineous Pakistani marriage; the parents of the cousin denied consanguinity. Agenesis of the corpus callosum, prominent root of the nose, large anteverted ears, facial hypertrichosis, small pupils with posterior synechiae, hypotonia, mild to moderate spastic palsy with hip dislocations, and hormonal dysfunction, presumably of hypothalamic origin, were other features. The children were almost blind, whether or not the cataracts had been operated on. The electroretinographic responses indicated dysfunction of both retinal rods and cones, and the visual evoked potentials confirmed optic nerve atrophy. The children were late walkers and were incontinent of urine and stools. In the differential diagnosis, Warburg et al. (1993) considered COFS syndrome (214150), CAMAK/CAMFAK syndromes (212540), Martsolf syndrome (212720), lethal Rutledge syndrome (270400), and lethal Neu-Laxova syndrome (256520).  ";
		//assertTrue(parser.tagInput(input).getMarkedText().contains("HP:0000007"));		
		
		input = "Argov and Yarom (1984) described a 'rimmed vacule myopathy' in Jews of Persian origin. The onset of this disorder usually occurred after the age of 20 years but before the middle of the fourth decade of life. Proximal and distal muscle weakness and wasting of the upper and lower limbs were progressive and resulted in severe incapacitation within 10 to 20 years. Despite this, there typically was sparing of the quadriceps muscles even in advanced stages of the disease, a feature unique to this form of inclusion body myopathy. It was not clear to the authors whether this disorder was primarily neurogenic or myopathic";
		assertTrue(parser.tagInput(input).getMarkedText().contains("HP:0003805"));
		
		input = "Osteogenesis imperfecta (see Byers, 1993) is characterized chiefly by multiple bone fractures, usually resulting from minimal trauma. Affected individuals have blue sclerae, normal teeth, and normal or near-normal stature (for growth curves, see Vetter et al., 1992). Fractures are rare in the neonatal period; fracture tendency is constant from childhood to puberty, decreases thereafter, and often increases following menopause in women and after the sixth decade in men. Fractures heal rapidly with evidence of a good callus formation, and, with good orthopedic care, without deformity. Hearing loss of conductive or mixed type occurs in about 50% of families, beginning in the late teens and leading, gradually, to profound deafness, tinnitus, and vertigo by the end of the fourth to fifth decade. Additional clinical findings may be thin, easily bruised skin, moderate joint hypermobility and kyphoscoliosis, hernias, and arcus senilis. Mitral valve prolapse, aortic valvular insufficiency, and a slightly larger than normal aortic root diameter have been identified in some individuals (Hortop et al., 1986), but it is not clear that these disorders are significantly more frequent than in the general population. Radiologically, wormian bones are common but bone morphology is generally normal at birth, although mild osteopenia and femoral bowing may be present. Vertebral body morphology in the adult is normal initially, but often develops the classic 'cod-fish' appearance (Steinmann et al., 1991).";
		assertTrue(parser.tagInput(input).getMarkedText().contains("HP:0000592"));
	}
	
}
