package genepi.check;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.StringTokenizer;

import genepi.io.table.TableReaderFactory;
import genepi.io.table.reader.CsvTableReader;
import genepi.io.table.reader.ITableReader;
import genepi.objects.CheckEntry;
import genepi.objects.HSDEntry;

public class HeteroplasmyBuilder {
	


	private String variantfile;
	private String outfile;
	private double vaf;


	public HeteroplasmyBuilder(String variantfile, String outfile) {
		this.variantfile = variantfile;
		this.outfile = outfile;
		this.vaf = 0.01;
	}

	public int build() throws MalformedURLException, IOException {
	
	FileInputStream fstream;
		try {
			
		//	input = directoryListing[0].getAbsolutePath();
			
			ITableReader idReader = TableReaderFactory.getReader(variantfile);


			HashMap<String, ArrayList<CheckEntry>> hm = new HashMap<String, ArrayList<CheckEntry>>();
			FileWriter fw;
			fw = new FileWriter(outfile);
			fw.write("SampleID\tRange\tHaplogroup\tPolymorphisms");
			fw.write(System.lineSeparator());
			try {
				while (idReader.next()) {
					CheckEntry entry = new CheckEntry();
					String id =  idReader.getString("SampleID"); //ID
					entry.setID(id);
					entry.setPOS(idReader.getInteger("Pos"));    //POS
					entry.setREF(idReader.getString("Ref"));	//REF
					entry.setALT(idReader.getString("Variant")); //ALT
					entry.setVAF(idReader.getDouble("Variant-Level")); //ALT

					if (hm.containsKey(id)) {
						hm.get(id).add(entry);
					} else if (hm.get(id) == null) {
						hm.put(id, new ArrayList<CheckEntry>());
						hm.get(id).add(entry);
					}
				}
				idReader.close();
			} catch (Exception e) {
				e.printStackTrace();
			}

			Iterator it = hm.entrySet().iterator();
			while (it.hasNext()) {
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


		} catch (Exception e) {
			System.out.println("ERROR");
			e.printStackTrace();
			return -1;
		}
		return 0;
	}
	
	
	
	

	public String getOutfile() {
		return outfile;
	}

	public void setOutfile(String outfile) {
		this.outfile = outfile;
	}

	public double getVaf() {
		return vaf;
	}

	public void setVaf(double vaf) {
		this.vaf = vaf;
	}
	
}
