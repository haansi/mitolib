package genepi.bam;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.util.StringTokenizer;

import genepi.base.Tool;


public class HaploCheckReader  extends Tool {

	public HaploCheckReader(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Split mitochondrial variants according heteroplasmy level - in profiles for HaploGrep 2\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input called variants with heteroplasmies");
		addParameter("out", "output file for HaploGrep 2");
		addOptionalParameter(
				"VAF",
				"optional: set the Variant Allele Frequency (VAF), default 0.1 = 10%)", DOUBLE);
	}

	@Override
	public int run() {

		String in = (String) getValue("in");
		String out = (String) getValue("out");
	
		double vaf = (Double) getValue("VAF");

		HaploCheckBuilder builder = new HaploCheckBuilder(in, out);
		builder.setVaf(vaf);


		try {
			return builder.build();
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			System.out.println("OOOJE");
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			
			System.out.println("OOOJE 2");
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
	}

	public static void main(String[] args) {
		new BAMReader(args).start();
	}

}
