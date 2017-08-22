package genepi.mitolib.haplogroup;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

import genepi.mitolib.objects.HeaderNames;

public class HaploMultiMap {

	Multimap<String, String> myMultimap = ArrayListMultimap.create();
	Multimap<String, String> phylotreemultimap = ArrayListMultimap.create();

	public void initMultiMaps() {
		// load data from phylotree - inverted
		readInverted(myMultimap);
		// load data from phylotre
		readPhylotree(phylotreemultimap);
	}

	public int getHaplogroup(String filename) {
		Long time = System.currentTimeMillis();

		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));

			String lineReader = br.readLine();
			if (lineReader.contains(HeaderNames.SampleId.colname()))
				lineReader = br.readLine(); // skip header
			while (lineReader != null) {

				String[] strsplit= lineReader.split("\t");
				
				String ID = strsplit[0]; // ID  1 Range 2 HG

				Collection<String> snpq = new ArrayList<String>();

				// add all haplogroups (values) to one SNP (key) from file starting at snp after tab 3
				for (int i =3; i< strsplit.length; i++) {
					String snp = strsplit[i];
					snpq.addAll(myMultimap.get(snp));
				}


				
				// generate a map with key = haplogroup and value size of query
				// SNPs found
				HashMap<String, Integer> mapq = new HashMap<String, Integer>();

				int max = 0;

				Iterator<String> iterator = snpq.iterator();
				while (iterator.hasNext()) {
					String hg = iterator.next();
					if (mapq.containsKey(hg)) {
						mapq.put(hg, mapq.get(hg) + 1);
						max = (max > (mapq.get(hg) + 1)) ? max : (mapq.get(hg));
					} else {
						mapq.put(hg, 1);
					}
				}

				Map.Entry<String, Integer> maxEntry = null;
				HashMap<String, Integer> mostfoundHGs = new HashMap<String, Integer>();

				// System.out.println("T check " + mostfoundHGs.contains("T"));

				for (Map.Entry<String, Integer> entry : mapq.entrySet()) {
					if (entry.getValue().compareTo(max) >= 0)
						mostfoundHGs.put(entry.getKey(), entry.getValue());
					if (maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) >= 0) {
						maxEntry = entry;
					}
				}

				// remove all entries that don't show the max. number of
				// possible SNPs
				if (maxEntry!=null){
				for (int i = 0; i < maxEntry.getValue(); i++) {
					mapq.values().removeAll(Collections.singleton(i));
				} 

				// no haplogroup with 1000 snps
				int min = 1000;
				String result = "";

				//find haplogroup with least differences 
				for (Map.Entry<String, Integer> entry : mapq.entrySet()) {
					String HG = entry.getKey();

					int sizeHG = phylotreemultimap.get(HG).size();
					if (sizeHG < min) {
						min = sizeHG;
						result = HG;
					}
				}
					System.out.println(ID + "\t" + result + "\t" + maxEntry.getValue() + "\t" + mapq.size());
				}
				else
				{
					System.out.println(ID + "\t" + "rCRS" + "\t" + 0 + "\t" + mapq.size());
				}
				
				lineReader = br.readLine();
			} // end while loop file entries
			br.close();

			Long time6 = System.currentTimeMillis();

			System.out.println("time " + (time6 - time) + " ms");

		} catch (Exception e) {
			e.printStackTrace();
			return -1;
		}

		return 0;
	}

	private static void readInverted(Multimap<String, String> myMultimap) {
		Long time = System.currentTimeMillis();

		String filename = "data/phylotree17.byPOS.txt";
		File file = new File(filename);
		BufferedReader br;

		String line = "-1";
		try {
			br = new BufferedReader(new FileReader(file));

			while ((line = br.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line, "\t,");
				String mtSNP = (String) st.nextElement();

				st.nextElement(); // skip Position
				st.nextElement(); // skip Info

				while (st.hasMoreTokens()) {
					String key = st.nextToken().trim();

					myMultimap.put(mtSNP, key);
				}
			}

		} catch (Exception e) {
			System.out.println("ERROR " + line);
			e.printStackTrace();
		}
		System.out.println("reading inverted phylotree in: " + (System.currentTimeMillis() - time) + " ms");
	}

	private static void readPhylotree(Multimap<String, String> phylotreemultimap) {
		Long time = System.currentTimeMillis();
		String filename = "data/phylotree17.hsd";
		File file = new File(filename);
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(file));

			String line;

			while ((line = br.readLine()) != null) {

				StringTokenizer st = new StringTokenizer(line, "\t,");
				String Haplogroup = (String) st.nextElement();
				st.nextElement(); // skip Range
				st.nextElement(); // skip Haplogroup

				while (st.hasMoreTokens()) {
					String key = st.nextToken().trim();

					phylotreemultimap.put(Haplogroup, key);
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("reading phylotree in: " + (System.currentTimeMillis() - time) + " ms");
	}

	public static void main(String[] args) {
		
		HaploMultiMap hm = new HaploMultiMap();
		hm.initMultiMaps();
		//hm.getHaplogroup("data/LG1_10_3000/LG1_10_3000.hsd");
		//hm.getHaplogroup("data/HG00740/HG00740.hsd");
		//hm.getHaplogroup("data/sim_H/sim_H2a2a1g.hsd");
		//hm.getHaplogroup("data/sim_H/sim_H.hsd");
		hm.getHaplogroup("data/1000G/1000G_P1.hsd");
		hm.getHaplogroup("data/1000G/1000G_P3.hsd");
		hm.getHaplogroup("data/LG1_10_3000/S4.hsd");
	}

}
