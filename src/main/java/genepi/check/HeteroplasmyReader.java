package genepi.check;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.util.StringTokenizer;

import genepi.base.Tool;


public class HeteroplasmyReader  extends Tool {

	public HeteroplasmyReader(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Split mitochondrial variants and heteroplasmies from mtDNA-Server (https://mtdna-server.uibk.ac.at) - in profiles for HaploGrep 2\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input Variants+Heteroplasmy from mtDNA-Server");
		addParameter("out", "output file for HaploGrep 2");
	}

	@Override
	public int run() {

		String in = (String) getValue("in");
		String out = (String) getValue("out");

		HeteroplasmyBuilder builder = new HeteroplasmyBuilder(in, out);

		try {
			return builder.build();
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
	}
}
