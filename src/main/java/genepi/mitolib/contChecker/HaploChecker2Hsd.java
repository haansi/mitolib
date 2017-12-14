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

import org.apache.commons.io.FilenameUtils;

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


public class HaploChecker2Hsd  extends Tool {

	public HaploChecker2Hsd(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Check mtDNA-Server raw data and generate input data for HaploGrep2 \n\n");

	}

	@Override
	public void createParameters() {
		addParameter("in", 	"input txt file");
		addParameter("out", "output folder");
		addParameter("VAF", "variant allele frequence", DOUBLE);
	}

	@Override
	public int run() {

		String input = (String) getValue("in");
		String out = (String) getValue("out");
		double inVaf = (Double) getValue("VAF");
		
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
			filename =  FilenameUtils.getBaseName(filename);
			
			createDirectory(out);
			
			CNVServer2Haplogrep splitter = new CNVServer2Haplogrep(emptyArgs);
			
			String infile = input;
			String fileHSD =out+File.separator+filename+".hsd";
			splitter.build(infile, fileHSD, vaf);
			
			System.out.println(fileHSD);

			System.out.println("Run-time: " + ((System.currentTimeMillis()-start)/1000)+ " sec");
			
		}catch (Exception e) {
					e.printStackTrace();
					return -1;
				}
		return 0;
		
	}
		
	
	private void createDirectory(String out) {
		new File(out).mkdirs();
		
	}


	public static void main(String[] args) {
		HaploChecker2Hsd h = new HaploChecker2Hsd(args);

		String input ="data/MixUps/rawM1-M4_0.001.txt";
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
