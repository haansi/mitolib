package mitolib;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.junit.Test;

import genepi.mitolib.splitter.HeteroplasmySplitter;


public class TestSplitter {

	@Test
	public void testSimulatedHGs() {
		
		String input = "data/sim_H.txt";
		String output = "data/sim_H.hsd";
		String[] args = new String[] {input, output};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(input, output, 0.01); //TEST WITH VAF 1% VAF
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		FileInputStream fis;
		try {
			fis = new FileInputStream(new File(output));
			String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
			fis.close();
			assertEquals("0a92b3289c05a4d2d14ee484e76f7c70", md5);
			} catch (FileNotFoundException e) {
			fail("file not generated");
			e.printStackTrace();
		} catch (IOException e) {
			fail("Hash not generated");
			e.printStackTrace();
		}
	
	}

}
