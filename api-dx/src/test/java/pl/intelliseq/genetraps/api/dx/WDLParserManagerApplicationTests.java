package pl.intelliseq.genetraps.api.dx;

import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.boot.test.autoconfigure.web.servlet.AutoConfigureMockMvc;
import org.springframework.boot.test.context.SpringBootTest;
import org.springframework.test.context.junit4.SpringRunner;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.request.MockMvcRequestBuilders;

import java.io.BufferedReader;
import java.io.FileReader;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.containsString;
import static pl.intelliseq.genetraps.api.dx.TestUser.PSYDUCK;

@RunWith(SpringRunner.class)
@AutoConfigureMockMvc
@SpringBootTest(webEnvironment = SpringBootTest.WebEnvironment.RANDOM_PORT)
public class WDLParserManagerApplicationTests {

	@Autowired
	private MockMvc mockMvc;

	@Test
	public void infoTest() {
		String response;
		BufferedReader reader;
		String testExpectedResult;
		try {
			reader = new BufferedReader(new FileReader("src/test/resources/basicTestResult"));
			testExpectedResult = reader.readLine();
			response = new ObjectMapper().readTree(
					mockMvc.perform(MockMvcRequestBuilders.get(String.format("/wdl/info"))
							.header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
							.andReturn().getResponse().getContentAsString()
			).toString();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		assertThat(response, containsString(testExpectedResult));
	}

	@Test
	public void infoSingleTest() {
		String response;
		BufferedReader reader;
		String testExpectedResult;
		try {
			reader = new BufferedReader(new FileReader("src/test/resources/testRes4"));
			testExpectedResult = reader.readLine();
			response = new ObjectMapper().readTree(
					mockMvc.perform(MockMvcRequestBuilders.get(String.format("/wdl/task-1/info"))
							.header("Authorization", "Bearer " + PSYDUCK.getAccessToken()))
							.andReturn().getResponse().getContentAsString()
			).toString();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		assertThat(response, containsString(testExpectedResult));
	}
}
