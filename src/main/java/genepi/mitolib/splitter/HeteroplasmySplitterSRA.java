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


public class HeteroplasmySplitterSRA  extends Tool {

	public HeteroplasmySplitterSRA(String[] args) {
		super(args);
	}
	@Override
	public void init() {

		System.out
				.println("Split mitochondrial pileup file from Peter SRA - in profiles for HaploGrep 2\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input raw file (pileup) from Peter");
		addParameter("vaf", "variant allele frequence ", DOUBLE);
	}
	
	
	@Override
	public int run() {

		String in = (String) getValue("in");
		double vaf = (double) getValue("vaf");

		try {
			return build(in, vaf);
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
	
	public int build(String variantfile, double vaf) throws MalformedURLException, IOException {
		try {
			
			//	input = directoryListing[0].getAbsolutePath();
				
				ITableReader idReader = TableReaderFactory.getReader(variantfile);


				HashMap<String, ArrayList<CheckEntry>> hm = new HashMap<String, ArrayList<CheckEntry>>();
				ArrayList<String> missingList = new ArrayList<>();
				int missings=0;
				int covAll=0;
				int linecount=0;
				String id =  variantfile; //SampleID
				
				
		try {
			while (idReader.next()) {
				CheckEntry entry = new CheckEntry();
				linecount++;
				
				entry.setID(id);
				
				entry.setPOS(idReader.getInteger("Position"));  //Pos
				String Ref = idReader.getString("Ref");
				entry.setREF(Ref);	//Ref
				
				int coverageUnfiltered= idReader.getInteger("Coverage");
				
				int RefCount= (int) (idReader.getInteger("RefCount"));
				
				int A= (int) (idReader.getInteger("A"));
				int C= (int) (idReader.getInteger("C"));
				int G= (int) (idReader.getInteger("G"));
				int T= (int) (idReader.getInteger("T"));
				
				int covTotal = A + C + G + T + RefCount;
				
				switch(Ref){
				case "A": A+=RefCount; break;
				case "C": C+=RefCount; break;
				case "G": G+=RefCount; break;
				case "T": T+=RefCount; break;
				}
				
				TreeMap<Integer, String> map = new TreeMap<Integer, String>(Collections.reverseOrder());
				
				map.put(A, "A");
				map.put(C, "C");
				map.put(G, "G");
				map.put(T, "T");

				if(covTotal>0){

			 covAll+=covTotal;
			 entry.setCOV(covTotal);
				String firstBase = (String) map.values().toArray()[0]; 
				double firstBasePerc=  Double.valueOf((map.keySet().toArray()[0])+"") / (covTotal); 
				
				String secondBase = (String) map.values().toArray()[1]; 
				double secondBasePerc= Double.valueOf( map.keySet().toArray()[1] +"") / (covTotal); 
		 
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
			    	entry.setVAF(firstBasePerc);
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
		 else{
			 missingList.add(entry.getPOS()+"");
			 missings++;
		 }
			}
			idReader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}	
		
		FileWriter fw;
		fw = new FileWriter(variantfile+".stats.txt");
	//	fw.write("SampleID\tlines\tQ30_filtered\tmeanCov_Q30\tsitesNotcovered");
		fw.write(System.lineSeparator());
		fw.write(id+ "\t"+linecount +"\t"+ missings  +"\t"+((covAll>0) ? (covAll/(linecount-missings)): 0)+"\t"+missingList.toString());
		fw.close();
		
		int counter = generateHSDfile(hm, variantfile, vaf);	
			
			System.out.println(counter +" samples processed\n" +counter*2 +" profiles written");

			} catch (Exception e) {
				System.out.println("ERROR");
				e.printStackTrace();
				return -1;
			}
			return 0;
		}

	public static int generateHSDfile(HashMap<String, ArrayList<CheckEntry>> hm, String outfile, double vaf) throws Exception{
		FileWriter fw, fwVariant;
		System.out.println("outfile " + outfile);
		fw = new FileWriter(outfile+".hsd");
		fwVariant = new FileWriter(outfile+".txt");
		
	//	fw.write("SampleID\tRange\tHaplogroup\tPolymorphisms");
		fw.write(System.lineSeparator());
		
		fwVariant.write("SampleID\tPos\tRef\tVariant\tVariant-Level\tCoverage-Total");
		fwVariant.write(System.lineSeparator());
		
		int counter =0;
		Iterator it = hm.entrySet().iterator();
		while (it.hasNext()) {
			counter++;
			Map.Entry pair = (Map.Entry) it.next();
			HSDEntry minor = new HSDEntry();
			HSDEntry major = new HSDEntry();
			minor.setID(pair.getKey() + "_maj");
			minor.setRANGE("1-16569");
			major.setID(pair.getKey() + "_min");
			major.setRANGE("1-16569");
			int hetcounter=0;
			ArrayList<CheckEntry> helpArray = hm.get(pair.getKey());
			for (int i = 0; i < helpArray.size(); i++) {
					if (helpArray.get(i).getREF().contains("-") || helpArray.get(i).getALT().contains("-")
						|| helpArray.get(i).getREF().equals("N") || helpArray.get(i).getALT().length() > 1
						|| helpArray.get(i).getREF().length() > 1) {
					// skip indel, and 3107 on rCRS;
				} else {
					if (helpArray.get(i).getVAF() < 0.5 ) {
						fwVariant.write(pair.getKey() + "\t" + helpArray.get(i).getPOS() + "\t" + helpArray.get(i).getREF() +  "\t" + helpArray.get(i).getALT() +  "\t" + helpArray.get(i).getVAF() +  "\t" + helpArray.get(i).getCOV() );
						fwVariant.write(System.lineSeparator());
						
						minor.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getREF());
						major.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getALT());
						hetcounter++;
					} else if (helpArray.get(i).getVAF() >= 0.5 && helpArray.get(i).getVAF() < 1-vaf){
						fwVariant.write(pair.getKey() + "\t" + helpArray.get(i).getPOS() + "\t" + helpArray.get(i).getREF() +  "\t" + helpArray.get(i).getALT() +  "\t" + helpArray.get(i).getVAF() +  "\t" + helpArray.get(i).getCOV());
						fwVariant.write(System.lineSeparator());
						
						minor.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getALT());
						major.appendPROFILES(helpArray.get(i).getPOS() + helpArray.get(i).getREF());
						hetcounter++;
					}
					else{
						fwVariant.write(pair.getKey() + "\t" + helpArray.get(i).getPOS() + "\t" + helpArray.get(i).getREF() +  "\t" + helpArray.get(i).getALT() +  "\t"+ helpArray.get(i).getVAF() +  "\t" + helpArray.get(i).getCOV());
						fwVariant.write(System.lineSeparator());
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
		fwVariant.close();
	return counter;
	}
	
	public static void main(String[] args){
		String in = "data/peter/test.mito";
		double vaf=0.01;
		
		HeteroplasmySplitterSRA hsra = new HeteroplasmySplitterSRA(args);
		try {
			hsra.build(in, vaf);
			in = "data/peter/ERR436061.mito";
			hsra.build(in, vaf);
			in = "data/peter/SRR1286931.mito";
			hsra.build(in, vaf);
			in = "data/peter/SRR1695629.mito";
			hsra.build(in, vaf);
			in = "data/peter/SRR924015.mito";
			hsra.build(in, vaf);
			in = "data/peter/SRR099317.mito";
			hsra.build(in, vaf);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
