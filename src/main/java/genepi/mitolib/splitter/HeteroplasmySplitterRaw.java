package genepi.mitolib.splitter;

import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import genepi.base.Tool;
import genepi.io.table.TableReaderFactory;
import genepi.io.table.reader.ITableReader;
import genepi.mitolib.objects.CheckEntry;
import genepi.mitolib.objects.HSDEntry;
import genepi.mitolib.objects.HeaderNames;


public class HeteroplasmySplitterRaw  extends Tool {

	public HeteroplasmySplitterRaw(String[] args) {
		super(args);
	}
	@Override
	public void init() {

		System.out
				.println("Split mitochondrial raw file from mtDNA-Server (https://mtdna-server.uibk.ac.at) - in profiles for HaploGrep 2\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input raw file (pileup) from mtDNA-Server");
		addParameter("vaf", "variant allele frequence ");
		addParameter("out", "output file for HaploGrep 2");
	}
	
	
	@Override
	public int run() {

		String in = (String) getValue("in");
		String out = (String) getValue("out");
		double vaf = 0.01;

		try {
			return build(in, out, vaf);
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
	
	public int build(String variantfile, String outfile, double vaf) throws MalformedURLException, IOException {
		try {
			
			//	input = directoryListing[0].getAbsolutePath();
				
				ITableReader idReader = TableReaderFactory.getReader(variantfile);


				HashMap<String, ArrayList<CheckEntry>> hm = new HashMap<String, ArrayList<CheckEntry>>();
	
		try {
			while (idReader.next()) {
				CheckEntry entry = new CheckEntry();
				String id =  idReader.getString("ID"); //SampleID
				entry.setID(id);
				
				entry.setPOS(idReader.getInteger("POS"));  //Pos
				entry.setREF(idReader.getString("rCRS"));	//Ref
				
				int covFWD= idReader.getInteger("COV-FWD");
				int covREV= idReader.getInteger("COV-REV");
				int A= (int) (idReader.getDouble("%A")*covFWD) +  (int)(idReader.getDouble("%a")*covREV);
				int C= (int) (idReader.getDouble("%C")*covFWD) +  (int)(idReader.getDouble("%c")*covREV);
				int G= (int) (idReader.getDouble("%G")*covFWD) +  (int)(idReader.getDouble("%g")*covREV);
				int T= (int) (idReader.getDouble("%T")*covFWD) +  (int)(idReader.getDouble("%t")*covREV);
				
				TreeMap<Integer, String> map = new TreeMap<Integer, String>(Collections.reverseOrder());
				
				map.put(A, "A");
				map.put(C, "C");
				map.put(G, "G");
				map.put(T, "T");
		 
				String firstBase = (String) map.values().toArray()[0]; 
				double firstBasePerc=  Double.valueOf((map.keySet().toArray()[0])+"") / (covFWD+covREV); 
				
				String secondBase = (String) map.values().toArray()[1]; 
				double secondBasePerc= Double.valueOf( map.keySet().toArray()[1] +"") / (covFWD+covREV); 
				
				boolean isVariant=false;
				if (firstBase.equals(entry.getREF()))
			    {
			    	entry.setALT(secondBase);
			    	entry.setVAF(secondBasePerc);
			    	if (secondBasePerc >= vaf){
			    		isVariant=true;
			    	}
			    		
			    }
			    else{
			    	entry.setALT(firstBase);
			    	entry.setVAF(1-firstBasePerc);
			    	isVariant=true;
			    }
		
				if (isVariant){
				if (hm.containsKey(id)) {
					hm.get(id).add(entry);
				} else if (hm.get(id) == null) {
					hm.put(id, new ArrayList<CheckEntry>());
					hm.get(id).add(entry);
				}
			}
			}
			idReader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}	
		
		
		
			int counter = generateHSDfile(hm, outfile, vaf);	
			
			System.out.println(counter +" samples processed\n" +counter*2 +" profiles written");

			} catch (Exception e) {
				System.out.println("ERROR");
				e.printStackTrace();
				return -1;
			}
			return 0;
		}

	public static int generateHSDfile(HashMap<String, ArrayList<CheckEntry>> hm, String outfile, double vaf) throws Exception{
		FileWriter fw;
		fw = new FileWriter(outfile);
		fw.write("SampleID\tRange\tHaplogroup\tPolymorphisms");
		fw.write(System.lineSeparator());
		int counter =0;
		Iterator it = hm.entrySet().iterator();
		while (it.hasNext()) {
			counter++;
			Map.Entry pair = (Map.Entry) it.next();
			HSDEntry minor = new HSDEntry();
			HSDEntry major = new HSDEntry();
			minor.setID(pair.getKey() + "_min");
			minor.setRANGE("1-16569");
			major.setID(pair.getKey() + "_maj");
			major.setRANGE("1-16569");
			int hetcounter=0;
			ArrayList<CheckEntry> helpArray = hm.get(pair.getKey());
			for (int i = 0; i < helpArray.size(); i++) {
				
				if (helpArray.get(i).getREF().contains("-") || helpArray.get(i).getALT().contains("-")
						|| helpArray.get(i).getREF().equals("N") || helpArray.get(i).getALT().length() > 1
						|| helpArray.get(i).getREF().length() > 1) {
					// skip indel, and 3107 on rCRS;
				} else {
					if (helpArray.get(i).getVAF() < 0.5) {
						minor.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getREF());
						major.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getALT());
						hetcounter++;
					} else if (helpArray.get(i).getVAF() >= 0.5 && helpArray.get(i).getVAF() < 1-vaf){
						minor.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getALT());
						major.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getREF());
						hetcounter++;
					}
					else{
						minor.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getALT());
						major.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getALT());
					}
					
				}
			}
			if (hetcounter>0){
			fw.write(minor.getString());
			fw.write(System.lineSeparator());
			fw.write(major.getString());
			fw.write(System.lineSeparator());
			}
			it.remove(); // avoids a ConcurrentModificationException
		}
		fw.close();
	return counter;
	}
}
