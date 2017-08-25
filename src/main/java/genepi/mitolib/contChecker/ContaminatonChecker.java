package genepi.mitolib.contChecker;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;

import genepi.base.Tool;
import genepi.io.table.TableReaderFactory;
import genepi.io.table.reader.CsvTableReader;
import genepi.io.table.reader.ITableReader;
import genepi.mitolib.objects.ContaminationEntry;
import genepi.mitolib.objects.HeaderNames;


public class ContaminatonChecker  extends Tool {

	public ContaminatonChecker(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Compare mitochondrial profiles from extended report in HaploGrep 2\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("inHG2", 	"input HaploGrep2 extended file");
		addParameter("inVar", 	"input variant file");
		addParameter("out", "output file of contaminated Samples");
	}

	@Override
	public int run() {

		String inHG2 = (String) getValue("inHG2");
		String inVar = (String) getValue("inVar");
		String out = (String) getValue("out");
	

		try {
			return build(inHG2, inVar, out);
		} catch (MalformedURLException e) {
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	
	public int build(String inHG2, String inVar, String outfile ) throws MalformedURLException, IOException {
		int countEntries=0;
		int countPossibleContaminated=0;
		int countContaminated=0;
		
		try {

			ITableReader readTableLevels = TableReaderFactory.getReader(inVar);
		
			//hahslevels contains key = sampleid+"-"+pos+variant. e.g.: HG00096-152C
			HashMap<String, Double> hashLevels = new HashMap<String, Double>();
			
			try {
				while (readTableLevels.next()) {
					hashLevels.put(readTableLevels.getString(HeaderNames.SampleId.colname())+"-"+readTableLevels.getString(HeaderNames.Position.colname())+readTableLevels.getString(HeaderNames.VariantBase.colname()), readTableLevels.getDouble(HeaderNames.VariantLevel.colname()));
				}
				readTableLevels.close();
			}catch (Exception e) {
				e.printStackTrace();
			}
			

			CsvTableReader readTableHaploGrep = new CsvTableReader(inHG2, '\t', true);
			NumberFormat formatter = new DecimalFormat("#0.00");  

			ArrayList<ContaminationEntry> contArray = new  ArrayList<ContaminationEntry>();
			FileWriter fw;
			fw = new FileWriter(new File(outfile));
			fw.write("SampleID\tContamination\tMinorHG\tMinorSNPs\tMinorLevel\tMajorHG\tMajorSNPs\tMajorLevel");
			fw.write(System.lineSeparator());
		
			try {
				while (readTableHaploGrep.next()) {
					ContaminationEntry centry = new ContaminationEntry();
					countEntries++;
					String id =  readTableHaploGrep.getString("SampleID"); //ID
					double weight =  readTableHaploGrep.getDouble("Overall_Rank"); //Rank
					centry.setSampleId(id.split("_maj")[0]);
					centry.setMajorId(readTableHaploGrep.getString("Haplogroup"));    //Major
			
					
					String notfound = readTableHaploGrep.getString("Not_Found_Polys");
					centry.setMajorRemaining(notfound.length() - notfound.replaceAll(" ", "").length());
					String majorfound = readTableHaploGrep.getString("Found_Polys");
					double meanMajor = getMeanScores(centry.getSampleId(), majorfound, hashLevels);
					
					//check second pair entry
			
					readTableHaploGrep.next();
					centry.setMinorId(readTableHaploGrep.getString("Haplogroup"));	  //Minor
					notfound = readTableHaploGrep.getString("Not_Found_Polys");
					centry.setMinorRemaining(notfound.length() - notfound.replaceAll(" ", "").length());
					String minorfound = readTableHaploGrep.getString("Found_Polys");
					double meanMinor = getMeanScores(centry.getSampleId(), minorfound, hashLevels);
					
					if (!centry.getMajorId().equals(centry.getMinorId())) {
						contArray.add(centry);
						countPossibleContaminated++;
						if (weight>0.7)
							{
							countContaminated ++;
							fw.write(centry.getSampleId()+"\tHigh\t"+centry.getMajorId() +"\t"+ formatter.format(meanMajor) +"\t"+ (majorfound.length() - majorfound.replaceAll(" ", "").length())+"\t" + centry.getMinorId()+"\t" +formatter.format(meanMinor) +"\t"+ (minorfound.length() - minorfound.replaceAll(" ", "").length())+"\n");
					}else if (notfound.length() - notfound.replaceAll(" ", "").length()>1){
							fw.write(centry.getSampleId()+"\tInc.\t"+centry.getMajorId() +"\t"+ formatter.format(meanMajor) +"\t"+ (majorfound.length() - majorfound.replaceAll(" ", "").length())+"\t" + centry.getMinorId()+"\t" +formatter.format(meanMinor) +"\t"+ (minorfound.length() - minorfound.replaceAll(" ", "").length())+"\n");
						}
					else	{
						fw.write(centry.getSampleId()+"\tPoss\t"+centry.getMajorId() +"\t"+ formatter.format(meanMajor) +"\t"+ (majorfound.length() - majorfound.replaceAll(" ", "").length())+"\t" + centry.getMinorId()+"\t"  +formatter.format(meanMinor) +"\t"+ (minorfound.length() - minorfound.replaceAll(" ", "").length())+"\n");
						
					}
					}
				
					else{
							
						//System.out.println(centry.getSampleId());
					}
				}
				readTableHaploGrep.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			System.out.println(countPossibleContaminated + " of " + countEntries + " Samples possibly contaminated" );
			System.out.println(countContaminated + " very likely" );
			
				
			fw.close();

		

		} catch (Exception e) {
			System.out.println("ERROR");
			e.printStackTrace();
			return -1;
		}
		// Everything fine
		return 0;
	}
	
	
	private double getMeanScores(String sampleId, String found, HashMap<String, Double> hmap) {
		
		double sum1 = 0;
		double sum2 =0;
		double stdev =0;

		int i=0;
		
		StringTokenizer st = new StringTokenizer(found, " ");
		while (st.hasMoreTokens()){
			String variant = st.nextToken();
			
			if (hmap.containsKey(sampleId+"-"+variant)){
				
				double value = hmap.get(sampleId+"-"+variant);
				if (value <0.99){
				sum1+= value;
				//sum2+=Math.pow(value, 2);
				//stdev = Math.sqrt(i*sum2 - Math.pow(sum1, 2))/i;
				//System.out.println(stdev);
				i++;
				}
			}
		}
		
	if (i>0){
		return sum1/i;
	}
	else
		return 0;
	}
	

	public static void main(String[] args) {
		new ContaminatonChecker(args).start();
	}

}
