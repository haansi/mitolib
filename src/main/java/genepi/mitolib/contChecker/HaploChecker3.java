package genepi.mitolib.contChecker;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.reflect.Array;
import java.net.MalformedURLException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import com.google.common.collect.Maps;

import genepi.base.Tool;
import genepi.io.table.TableReaderFactory;
import genepi.io.table.reader.CsvTableReader;
import genepi.io.table.reader.ITableReader;
import genepi.mitolib.bam.BAMReader;
import genepi.mitolib.haplogroup.HaploGrepCMD;
import genepi.mitolib.objects.ContaminationEntry;
import genepi.mitolib.objects.HeaderNames;
import genepi.mitolib.splitter.CNVServer2Haplogrep;
import genepi.mitolib.splitter.HeteroplasmySplitter;
import server.main.HaplogrepCMD;


public class HaploChecker3  extends Tool {

	public HaploChecker3(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Check mtDNA-Server variant data for contamination detection based on mitochondrial phylogeny (Phylotree 17)\n\n");

	}

	@Override
	public void createParameters() {
		addParameter("in", 	"input txt file");
		addParameter("vaf", "variant allele frequence", DOUBLE);
		addParameter("out", "output folder");
	}

	@Override
	public int run() {

		String input = (String) getValue("in");
		double inVaf = (Double) getValue("vaf");
		String out = (String) getValue("out");
	
		try {
			return build(input, inVaf, out);
		} catch (MalformedURLException e) {
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	
	
	public int build(String input, double vaf, String out) throws MalformedURLException, IOException {
		String[] emptyArgs=new String[]{""};
		try {
			double start = System.currentTimeMillis();
			
			String filename= input.split(File.separator)[input.split(File.separator).length-1];
			
			HeteroplasmySplitter splitter = new HeteroplasmySplitter(emptyArgs);
			
			String infile = input;
			String fileHSD =out+File.separator+filename+".hsd";
			//String fileVariants =out+File.separator+filename+".txt";
			System.out.println("HE$RE we are " + infile + " " +fileHSD + " " + vaf);
			
			String fileHaploGrep =out+File.separator+filename+".haplogrep.txt";
			splitter.build(infile, fileHSD, vaf);
		
			String[] testArgs ={"-format", "hsd", "-in", fileHSD,"-out",fileHaploGrep, "-phylotree","17"};
			
			HaploGrepCMD hg = new HaploGrepCMD(testArgs);
			hg.setArgs(testArgs);
			hg.run();
		
			
			String fileContamination =out+File.separator+filename+".contamination.txt";
			ContaminatonChecker contChecker = new ContaminatonChecker(emptyArgs);
			contChecker.build(fileHaploGrep, infile, fileContamination, vaf, null);
		
			System.out.println("Run-time: " + ((System.currentTimeMillis()-start)/1000)+ " sec");
			
		}catch (Exception e) {
					e.printStackTrace();
					return -1;
				}
		return 0;
		
	}
		
	
	public static void main(String[] args) {
		HaploChecker3 h = new HaploChecker3(args);

		String input ="data/MixUps/variantsLocal.txt";
		String out="data/output";
		double vaf =0.01;
		try {
			h.build(input, vaf, out);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	System.exit(0);	
	}

}
