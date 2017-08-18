package mitolib;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.junit.Test;

import genepi.mitolib.contChecker.ContaminatonChecker;

public class TestContChecker {

	@Test
	public void testSimulated_H() {
		String inputHaploGrep = "data/sim_H_haplogrep.txt";
		String inputLevels = "data/sim_H.txt";
		String output = "data/sim_H_contaminated.txt";
		String[] args = new String[] {inputHaploGrep, inputLevels, output};
		ContaminatonChecker contChecker = new ContaminatonChecker(args);
		try {
			contChecker.build(inputHaploGrep, inputLevels, output);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		FileInputStream fis;
		try {
			fis = new FileInputStream(new File(output));
			String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
			fis.close();
			assertEquals("861f1cfab86f2ed49f69bf364f3c3d21", md5);
			} catch (FileNotFoundException e) {
			fail("file not generated");
			e.printStackTrace();
		} catch (IOException e) {
			fail("Hash not generated");
			e.printStackTrace();
		}
	
	}

}
