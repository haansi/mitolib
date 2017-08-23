package mitolib;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import genepi.mitolib.contChecker.ContaminatonChecker;
import genepi.mitolib.haplogroup.HaploGrepCMD;
import genepi.mitolib.splitter.HeteroplasmySplitter;
import server.main.HaplogrepCMD;

public class TestContChecker {

	@Test
	public void testSimulated_H() {
		
		String inputLevels = "data/sim_H/sim_H.txt";
		String outputHsd = "data/sim_H/sim_H.hsd";
		String outputHaploGrep = "data/sim_H/sim_H_haplogrep.txt";
		String output = "data/sim_H/sim_H_contaminated.txt";

	
		String[] args = new String[] {inputLevels, outputHsd};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(inputLevels, outputHsd, 0.01); //TEST WITH VAF 1% VAF
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		String[] testArgs ={"-format", "hsd", "-in", outputHsd,"-out",outputHaploGrep, "-phylotree","17"};
		String strArgs ="--format hsd --in "+ outputHsd + "--out "+ outputHaploGrep+ "--phylotree 17";
		
		HaploGrepCMD hg = new HaploGrepCMD(testArgs);
		hg.setArgs(testArgs);
		hg.run();
		
		args = new String[] {outputHaploGrep, inputLevels, output};

		ContaminatonChecker contChecker = new ContaminatonChecker(args);
		try {
			contChecker.build(outputHaploGrep, inputLevels, output);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		FileInputStream fis;
		try {
			fis = new FileInputStream(new File(output));
			String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
			fis.close();
			assertEquals("ab128f17d5f4a0db97abd2c2a35196e9", md5);
			} catch (FileNotFoundException e) {
			fail("file not generated");
			e.printStackTrace();
		} catch (IOException e) {
			fail("Hash not generated");
			e.printStackTrace();
		}
	
	}
	
	
	
	@Test
	public void test_HG00740() {
		
		String inputLevels = "data/HG00740/HG00740.txt";
		String outputHsd = "data/HG00740/HG00740.hsd";
		String outputHaploGrep = "data/HG00740/HG00740_haplogrep.txt";
		String output = "data/HG00740/HG00740_contaminated.txt";

	
		String[] args = new String[] {inputLevels, outputHsd};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(inputLevels, outputHsd, 0.01); //TEST WITH VAF 1% VAF
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		String[] testArgs ={"-format", "hsd", "-in", outputHsd,"-out",outputHaploGrep, "-phylotree","17"};
		String strArgs ="--format hsd --in "+ outputHsd + "--out "+ outputHaploGrep+ "--phylotree 17";
		
		HaploGrepCMD hg = new HaploGrepCMD(testArgs);
		hg.setArgs(testArgs);
		hg.run();
		
		args = new String[] {outputHaploGrep, inputLevels, output};

		ContaminatonChecker contChecker = new ContaminatonChecker(args);
		try {
			contChecker.build(outputHaploGrep, inputLevels, output);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		FileInputStream fis;
		try {
			fis = new FileInputStream(new File(output));
			String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
			fis.close();

			assertEquals("00975ba72190371262d96a9a59516bf1", md5);
			} catch (FileNotFoundException e) {
			fail("file not generated");
			e.printStackTrace();
		} catch (IOException e) {
			fail("Hash not generated");
			e.printStackTrace();
		}
	
	}

	@Test
	public void testLG1_10_3000() {
		System.out.println("LG1-10");
		String inputLevels = "data/LG1_10_3000/LG1_10_3000.txt";
		String outputHsd = "data/LG1_10_3000/LG1_10_3000.hsd";
		String outputHaploGrep = "data/LG1_10_3000/LG1_10_3000_haplogrep.txt";
		String output = "data/LG1_10_3000/LG1_10_3000_contaminated.txt";

	
		String[] args = new String[] {inputLevels, outputHsd};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(inputLevels, outputHsd, 0.01); //TEST WITH VAF 1% VAF
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		String[] testArgs ={"-format", "hsd", "-in", outputHsd,"-out",outputHaploGrep, "-phylotree","17"};
		String strArgs ="--format hsd --in "+ outputHsd + "--out "+ outputHaploGrep+ "--phylotree 17";
		
		HaploGrepCMD hg = new HaploGrepCMD(testArgs);
		hg.setArgs(testArgs);
		hg.run();
		
		args = new String[] {outputHaploGrep, inputLevels, output};

		ContaminatonChecker contChecker = new ContaminatonChecker(args);
		try {
			contChecker.build(outputHaploGrep, inputLevels, output);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		FileInputStream fis;
		try {
			fis = new FileInputStream(new File(output));
			String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
			fis.close();

			assertEquals("00975ba72190371262d96a9a59516bf1", md5);
			} catch (FileNotFoundException e) {
			fail("file not generated");
			e.printStackTrace();
		} catch (IOException e) {
			fail("Hash not generated");
			e.printStackTrace();
		}
	
	}
	
	
	@Test
	public void testLG1_10() {
		System.out.println("LG1-10 ");
		String inputLevels = "/home/hansi/git/greenVC/output/LG-1-10-NEB.bam.txt";
		String outputHsd = "/home/hansi/git/greenVC/output/LG-1-10-NEB.hsd";
		String outputHaploGrep = "/home/hansi/git/greenVC/output/LG-1-10-haplogrep.txt";
		String output = "/home/hansi/git/greenVC/output/LG-1-10-NEB.contaminated.txt";

	
		String[] args = new String[] {inputLevels, outputHsd};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(inputLevels, outputHsd, 0.01); //TEST WITH VAF 1% VAF
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		String[] testArgs ={"-format", "hsd", "-in", outputHsd,"-out",outputHaploGrep, "-phylotree","17"};
		String strArgs ="--format hsd --in "+ outputHsd + "--out "+ outputHaploGrep+ "--phylotree 17";
		
		HaploGrepCMD hg = new HaploGrepCMD(testArgs);
		hg.setArgs(testArgs);
		hg.run();
		
		args = new String[] {outputHaploGrep, inputLevels, output};

		ContaminatonChecker contChecker = new ContaminatonChecker(args);
		try {
			contChecker.build(outputHaploGrep, inputLevels, output);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		FileInputStream fis;
		try {
			fis = new FileInputStream(new File(output));
			String md5 = org.apache.commons.codec.digest.DigestUtils.md5Hex(fis);
			fis.close();

			assertEquals("00975ba72190371262d96a9a59516bf1", md5);
			} catch (FileNotFoundException e) {
			fail("file not generated");
			e.printStackTrace();
		} catch (IOException e) {
			fail("Hash not generated");
			e.printStackTrace();
		}
	
	}
	
	
	
	@Test
	public void test_S4() {
		System.out.println("S4");

		String outputHsd = "data/LG1_10_3000/seb.hsd";
		String outputHaploGrep = "data/LG1_10_3000/seb.txt";

		
		String[] testArgs ={"-format", "hsd", "-in", outputHsd,"-out",outputHaploGrep, "-phylotree","17"};
		HaploGrepCMD hg = new HaploGrepCMD(testArgs);
		
		hg.setArgs(testArgs);
		hg.run();

	
	}
	
}
