package mitolib;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import org.junit.Test;

import genepi.io.table.TableReaderFactory;
import genepi.io.table.reader.ITableReader;
import genepi.mitolib.contChecker.ContaminatonChecker;
import genepi.mitolib.haplogroup.HaploGrepCMD;
import genepi.mitolib.objects.HeaderNames;
import genepi.mitolib.splitter.HeteroplasmySplitter;


public class TestContChecker {
	double vaf =0.01;
	@Test
	public void testSimulated_H() {
		
		String inputLevels = "data/sim_H/sim_H.txt";
		String outputHsd = "data/sim_H/sim_H.hsd";
		String outputHaploGrep = "data/sim_H/sim_H_haplogrep.txt";
		String output = "data/sim_H/sim_H_contaminated.txt";
		double vaf =0.01;
	
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
			contChecker.build(outputHaploGrep, inputLevels, output,vaf, null);
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
	public void test1000G_greenVC() {
		
		String path="data/1000G/1000G_txt/";
		String inputLevels = path+"all.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
		
	
		String[] args = new String[] {inputLevels, outputHsd};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(inputLevels, outputHsd, vaf); //TEST WITH VAF 1% VAF
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000G_High_chip_Lofreq() {
		
		String path="data/validatebamid//";
		String inputLevels = path+"High_Chip_Mix_Lofreq.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
		
	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000G_High_free_mtDNA_server() {
		
		String path="data/validatebamid//";
		String inputLevels = path+"High_Free_Mix_mtdna_server.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";
		double vaf =0.01;

		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
		
	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000G_High_chip_mix_haplocheck_server() {
		
		String path="data/validatebamid//";
		String inputLevels = path+"High_Chip_Mix_cnv_server.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
		
	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000G_High_free_Lofreq() {
		
		String path="data/validatebamid//";
		String inputLevels = path+"High_Free_Mix_Lofreq.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
		
	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000GP3_0_5() {
		
	
		String path="data/1000G/";
		String inputLevels = path+"1000GP3_mtDNA_server_0_5_noBAQ.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+".contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
	
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
			contChecker.build(outputHaploGrep, inputLevels, output,vaf, verifyBamIDScoresFree);
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
	public void test1000GP3_H() {
		
	
		String path="data/1000G/";
		String inputLevels = path+"1000Gp3.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000GP3_BAQ() {
		
	
		String path="data/1000G/BAQ/";
		String inputLevels = path+"variantsLocal.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+".contaminated.txt";
		double vaf =0.01;
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont_max.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("chip_contam");
				verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, verifyBamIDScoresFree);
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
	public void test1000_HG00s() {
		
	
		String path="data/1000G/";
		String inputLevels = path+"HG00119_HG00739_HG00740.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contaminated.txt";


	
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
		long time = System.currentTimeMillis();
		ContaminatonChecker contChecker = new ContaminatonChecker(args);
		try {
			contChecker.build(outputHaploGrep, inputLevels, output,vaf, null);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
System.out.println("Run time: " +(System.currentTimeMillis() -time) );
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, null);
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf, null);
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
	public void testVerifyBamID_greenVC0_1_8() {
		System.out.println("VerifyBamIDs greenVC");
		
		String path="data/validatebamid/";
		String inputLevels = path+"1000G_greenVC.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contamintated.txt";
		double vaf =0.007; 
		
		HashMap<String, Double> verifyBamIDScoresFree = new HashMap<String, Double>();
		HashMap<String, Double> verifyBamIDScoresChip = new HashMap<String, Double>();
		
		ITableReader readTableLevels = TableReaderFactory.getReader("data/1000G/VerifyBamId_cont.txt");
		try {
			while (readTableLevels.next()) {
				String SampleID = readTableLevels.getString("ID");
				//double free_contam = readTableLevels.getDouble("free_contam");
				double chip_contam = readTableLevels.getDouble("contam");
				//verifyBamIDScoresFree.put(SampleID, free_contam);
				verifyBamIDScoresChip.put(SampleID, chip_contam);
			}
		}
			catch (Exception e) {
				// TODO: handle exception
			}
		
	
		String[] args = new String[] {inputLevels, outputHsd};
		HeteroplasmySplitter splitter = new HeteroplasmySplitter(args);
		try {
			splitter.build(inputLevels, outputHsd, vaf); //TEST WITH VAF 1% VAF
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
			contChecker.build(outputHaploGrep, inputLevels, output,vaf, verifyBamIDScoresChip);
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
	public void testVerifyBamID_mtDNA_server0_9() {
		System.out.println("VerifyBamIDs mtdna_server");
		
		String path="data/validatebamid/";
		String inputLevels = path+"1000G_mtDNA-Server0_9.txt";
		String outputHsd = inputLevels.split("\\.")[0]+".hsd";
		String outputHaploGrep = inputLevels.split("\\.")[0]+".haplogrep.txt";
		String output = inputLevels.split("\\.")[0]+"contamintated.txt";

	
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf,null);
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
			contChecker.build(outputHaploGrep, inputLevels, output, vaf,null);
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
