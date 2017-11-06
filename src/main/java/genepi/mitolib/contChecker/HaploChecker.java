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
import genepi.mitolib.splitter.HeteroplasmySplitter;
import server.main.HaplogrepCMD;


public class HaploChecker  extends Tool {

	public HaploChecker(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Checkes if a BAM file shows signs of contamination based on the mitochondrial phylogeny (Phylotree 17)\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input BAM file");
		addParameter("ref", "reference FASTA file");
		
		addParameter("out", "output folder");
	}

	@Override
	public int run() {

		String inBam = (String) getValue("in");
		String inRef = (String) getValue("ref");
		String out = (String) getValue("out");
	
		try {
			return build(inBam, inRef, out);
		} catch (MalformedURLException e) {
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	
	
	public int build(String inBam, String ref, String out) throws MalformedURLException, IOException {
		double vaf =0.01;
		int mapqual =200;
		int qual =20;
		String[] emptyArgs=new String[]{""};
		try {
			double start = System.currentTimeMillis();
			BAMReader br = new BAMReader(new String[]{""});
			br.build(inBam, out, ref, vaf, qual, mapqual);
			
			
			HeteroplasmySplitter splitter = new HeteroplasmySplitter(emptyArgs);
			
			String[] tokens = inBam.split(File.separator);
			String filenam = tokens[tokens.length-1];
			
			String infile = out+File.separator+filenam+".txt";
			String fileHSD =out+File.separator+filenam+".hsd";
			
			String fileHaploGrep =out+File.separator+filenam+".haplogrep.txt";
			splitter.build(infile, fileHSD, vaf);
			
			String[] testArgs ={"-format", "hsd", "-in", fileHSD,"-out",fileHaploGrep, "-phylotree","17"};
			
			HaploGrepCMD hg = new HaploGrepCMD(testArgs);
			hg.setArgs(testArgs);
			hg.run();
		
			
			String fileContamination =out+File.separator+filenam+".contamination.txt";
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
		HaploChecker h = new HaploChecker(args);
		//String inBam = "data/HG01500.IBS.exome.MT.bam";
		String inBam ="data/HG00740/HG00740.mapped.ILLUMINA.bwa.PUR.low_coverage.20120522.bam";
		String out="data/output";
		String ref ="data/rcrs.fasta";
		try {
			h.build(inBam, ref, out);
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
