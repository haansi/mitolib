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
			return build(inHG2, inVar, out, null);
		} catch (MalformedURLException e) {
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			e.printStackTrace();
			return -1;
		}
	}
	
	
	
	public int build(String inHG2, String inVar, String outfile, HashMap<String, Double> verifyBam) throws MalformedURLException, IOException {
		double countEntries=0;
		double countPossibleContaminated=0;
		double countContaminated=0;
		
		try {

			ITableReader readTableLevels = TableReaderFactory.getReader(inVar);
		
			HashMap<String,ArrayList<Integer>> coverageMap=new HashMap<String,ArrayList<Integer>>(); 
			//hahslevels contains key = sampleid+"-"+pos+variant. e.g.: HG00096-152C
			HashMap<String, Double> heteroLevels = new HashMap<String, Double>();
			HashMap<String, Integer> homoplasmies = new HashMap<String, Integer>();
			HashMap<String, Integer> homoplasmiesMeta = new HashMap<String, Integer>();
			try {
				while (readTableLevels.next()) {
					double vaf = readTableLevels.getDouble(HeaderNames.VariantLevel.colname());
					String ID = readTableLevels.getString(HeaderNames.SampleId.colname());
					String key = ID+"-"+readTableLevels.getString(HeaderNames.Position.colname())+readTableLevels.getString(HeaderNames.VariantBase.colname());
					double value = readTableLevels.getDouble(HeaderNames.VariantLevel.colname());
					int cov= readTableLevels.getInteger("Coverage-Total");
					if (vaf<1-0.01)
					{
						heteroLevels.put(key,value);
					}
					else{
						if (homoplasmiesMeta.containsKey(ID))
							homoplasmiesMeta.put(ID, homoplasmiesMeta.get(ID)+1);
						else
							homoplasmiesMeta.put(ID, 1);
						
						homoplasmies.put(key, 1);
					}
					if (coverageMap.get(ID) == null) {
						coverageMap.put(ID, new ArrayList<Integer>());
					}
					coverageMap.get(ID).add(cov);
				}
				readTableLevels.close();
			}catch (Exception e) {
				e.printStackTrace();
			}
			

			CsvTableReader readTableHaploGrep = new CsvTableReader(inHG2, '\t', true);
			NumberFormat formatter = new DecimalFormat("#0.00");  

			ArrayList<ContaminationEntry> contArray = new  ArrayList<ContaminationEntry>();
			FileWriter fw = new FileWriter(new File(outfile));
			fw.write("SampleID\tContamination\tMinorHG\tMinorSNPs\tMinorLevel\tMinorHGvariants\tMajorHG\tMajorSNPs\tMajorLevel\tMajorHGvariants\tVerifyScore\tmeanCovVar");
			fw.write(System.lineSeparator());
		
			try {
				while (readTableHaploGrep.next()) {
					ContaminationEntry centry = new ContaminationEntry();
					countEntries++;
					String id =  readTableHaploGrep.getString("SampleID"); //ID
					double verifyScore=0;
					try{
					if (verifyBam.containsKey(id.substring(0,7)));
						verifyScore=verifyBam.get(id.substring(0,7));
						
					} catch (Exception e) {
						// TODO: handle exception
					}
						
					double weight =  readTableHaploGrep.getDouble("Overall_Rank"); //Rank
					centry.setSampleId(id.split("_maj")[0]);
					centry.setMajorId(readTableHaploGrep.getString("Haplogroup"));    //Major
			
					
					String notfound = readTableHaploGrep.getString("Not_Found_Polys");
					centry.setMajorRemaining(notfound.length() - notfound.replaceAll(" ", "").length());
					String majorfound = readTableHaploGrep.getString("Found_Polys");
					double meanMajor = getMeanScores(centry.getSampleId(), majorfound, heteroLevels);
					int[] countHomoplMajor=countHomoplasmies(centry.getSampleId(), majorfound, homoplasmies, homoplasmiesMeta);
					double meanCov = getMeaCoverage(id.split("_maj")[0], coverageMap);
					//check second pair entry
			
					readTableHaploGrep.next();
					centry.setMinorId(readTableHaploGrep.getString("Haplogroup"));	  //Minor
					notfound = readTableHaploGrep.getString("Not_Found_Polys");
					centry.setMinorRemaining(notfound.length() - notfound.replaceAll(" ", "").length());
					String minorfound = readTableHaploGrep.getString("Found_Polys");
					double meanMinor = getMeanScores(centry.getSampleId(), minorfound, heteroLevels);
					int[] countHomoplMinor=countHomoplasmies(centry.getSampleId(), minorfound, homoplasmies, homoplasmiesMeta);
					
					int majMutfound = majorfound.length() - majorfound.replaceAll(" ", "").length();
					int minMutfound = minorfound.length() - minorfound.replaceAll(" ", "").length();
					
					String homoplMajor = countHomoplMajor[0]+"/"+countHomoplMajor[1];
					String homoplMinor = countHomoplMinor[0] +"/"+countHomoplMinor[1];
				
					//check if Haplogroup names are different:
					if (!centry.getMajorId().equals(centry.getMinorId())) {
						contArray.add(centry);
						countPossibleContaminated++;
							
						//check if one of the haplogroups is defined by at least 2 heteroplasmic variants
						if ((majMutfound - countHomoplMajor[0]) > 2 ||  (minMutfound - countHomoplMinor[0]) >2) {
							countContaminated++;
							fw.write(centry.getSampleId() + "\tHigh\t" + centry.getMajorId() + "\t"
									+ formatter.format(meanMajor) + "\t" + homoplMajor + "\t"
									+ (majMutfound - countHomoplMajor[0]) + "\t" + centry.getMinorId() + "\t"
									+ formatter.format(meanMinor) + "\t" + homoplMinor + "\t"
									+ (minMutfound - countHomoplMinor[0]) +"\t"+verifyScore + "\t"+meanCov+"\n");
						} else if ((minMutfound - countHomoplMinor[0]) >1) {// (notfound.length()
																				// -
																				// notfound.replaceAll("
																				// ",
																				// "").length()>1){
							fw.write(centry.getSampleId() + "\tPoss\t" + centry.getMajorId() + "\t"									+ formatter.format(meanMajor) + "\t" + homoplMajor + "\t"
									+ (majMutfound - countHomoplMajor[0]) + "\t" + centry.getMinorId() + "\t"
									+ formatter.format(meanMinor) + "\t" + homoplMinor + "\t"
									+ (minMutfound - countHomoplMinor[0]) +"\t"+verifyScore + "\t"+meanCov+ "\n");
						} else {
							fw.write(centry.getSampleId() + "\tPoss\t" + centry.getMajorId() + "\t"
									+ formatter.format(meanMajor) + "\t" + homoplMajor + "\t"
									+ (majMutfound - countHomoplMajor[0]) + "\t" + centry.getMinorId() + "\t"
									+ formatter.format(meanMinor) + "\t" + homoplMinor + "\t"
									+ (minMutfound - countHomoplMinor[0]) +"\t"+verifyScore + "\t"+meanCov+ "\n");
						}

					}

					else {
						fw.write(centry.getSampleId() + "\tNone\t" + centry.getMajorId() + "\t"
								+ formatter.format(meanMajor) + "\t" + homoplMajor + "\t"
								+ (majMutfound - countHomoplMajor[0]) + "\t" + centry.getMinorId() + "\t"
								+ formatter.format(meanMinor) + "\t" + homoplMinor + "\t"
								+ (minMutfound - countHomoplMinor[0]) +"\t"+verifyScore +"\t"+meanCov+ "\n");
						
						// samples where no contamination was found
						// System.out.println(centry.getSampleId());
					}
				}
				readTableHaploGrep.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			
			FileWriter fwMeta = new FileWriter(new File(outfile+".meta"));
			fwMeta.write("group\tvalue\tnumbers\n");
			fwMeta.write("High\t"+countContaminated/countEntries + "\t" + countContaminated+"\n");
			fwMeta.write("Possible\t"+(countPossibleContaminated-countContaminated)/countEntries + "\t" + (countPossibleContaminated-countContaminated)+"\n");
			fwMeta.write("Not found\t"+(countEntries- countPossibleContaminated)/countEntries + "\t" + (countEntries- countPossibleContaminated)+"\n");
			fwMeta.close();
			
			
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
				/*sum2+=Math.pow(value, 2);
				stdev = Math.sqrt(i*sum2 - Math.pow(sum1, 2))/i;
				System.out.println(stdev); */
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
	

	
	private int[] countHomoplasmies(String sampleId, String found, HashMap<String, Integer> hmap, HashMap<String, Integer> hmapSize) {

		int[] result = new int[2]; //0 = homoplasmies in haplogroup found 
								   //1 = all homoplasmies in this sample	
	
			
		HashMap<String, Integer> helpMap= new HashMap<>();

		StringTokenizer st = new StringTokenizer(found, " ");
		int common=0;
		while (st.hasMoreTokens()){
			String variant = st.nextToken();
			String key= sampleId+"-"+variant;
			helpMap.put(key, 1);
			if (hmap.containsKey(key))
				common++;
		}
		long start=System.currentTimeMillis();
		 result[0]= common; //Maps.difference(hmap, helpMap).entriesInCommon().size();
		// System.out.println(System.currentTimeMillis()-start);
		 result[1]=hmapSize.get(sampleId);
		return result;
	}
	
	
	private double getMeaCoverage(String sampleId, HashMap<String, ArrayList<Integer>> covMap) {
		ArrayList<Integer> entries = covMap.get(sampleId);
		int sum=0;
		
		for (int i=0; i < entries.size(); i++){
			sum+=entries.get(i);
		}
		
		return sum/entries.size();
	}
	
		
	
	public static void main(String[] args) {
		new ContaminatonChecker(args).start();
	}

}
