package genepi.mitolib.vcf;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import com.google.common.io.FileWriteMode;

import genepi.base.Tool;
import genepi.mitolib.objects.CheckEntry;
import genepi.mitolib.objects.HeaderNames;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFContigHeaderLine;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;


public class VCFReader extends Tool {
	


	public VCFReader(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Split mitochondrial variants according the VCF file generated with MToolBox\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input called variants with MToolBox");
		addParameter("out", "output files for HaploGrep 2");
	}

	@Override
	public int run() {

		String in = (String) getValue("in");
		String out = (String) getValue("out");
	
		try {
			return build(in, out);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			System.out.println("Something somewhere went terribly wrong");
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			
			System.out.println("Something somewhere went terribly wrong, again");
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
	}

	public int build(String vcffile, String outfile) throws MalformedURLException, IOException {
		ArrayList<String> columnIDs = new ArrayList<>();
			try {
				boolean requireIndex = false;
				final VCFFileReader reader = new VCFFileReader(new File(vcffile), requireIndex);
				FileWriter fw = new FileWriter(new File(outfile));
				ArrayList<String> input = reader.getFileHeader()
						.getSampleNamesInOrder();
				Iterator<VCFContigHeaderLine> headerContigIt = reader.getFileHeader().getContigLines().iterator(); 
				Iterator<VCFFormatHeaderLine> headerFormatIt = reader.getFileHeader().getFormatHeaderLines().iterator();

				String ID ="";
				String FORMAT ="";
				
				
				CloseableIterator<VariantContext> it = reader.iterator();
				ArrayList<String> result = new ArrayList<String>();
				// write sample IDS + basic HSD format
				for (int i = 0; i < input.size(); i++) {
					columnIDs.add(input.get(i));
				}

				while (it.hasNext()) {
					/* get next variation and save it */
					VariantContext vc = it.next();
				
					for (int i = 0; i < input.size(); i++) {
						String genotype = vc.getGenotype(input.get(i))
								.getGenotypeString();
					

						if (genotype.contains("/")) // TODO later check for heteroplasmy
							genotype = genotype.split("/")[0];
						
				
						
						String type = vc.getGenotype(input.get(i)).getType().toString();
						String ref = vc.getReference().toString().replace("*", "");

						if (vc.getGenotype(input.get(i)).isHet())
							{
						
							String HF = vc.getGenotype(input.get(i)).getAnyAttribute("HF")+"";
							if (HF.contains(","))
							{
								result.add(HF.split(",")[(HF.length() - HF.replace(",","").length())]);
							}
							
						else
						{
							 //System.out.println(input.get(i) + "\t" +vc.getStart() + "  "+ "" + " " + genotype + " " + type + " " + ref);
							result.add(HF);

						}
							}
							else{
								result.add(0+"");
							}
						
				}
				}
	
				for (int i =0; i < result.size(); i++){
					for (int j =0; j < columnIDs.size(); j++)
						fw.append(result.get(i)+"\t");
				fw.append("\n");
				}
				
			reader.close();
			it.close();
				
			fw.close();
				
		
				System.out.println("File written");
			} catch (Exception e) {
				System.out.println("ERROR");
				e.printStackTrace();
			} 
			// Everything fine
			return 0;
		}
	
	public static void main(String args[]){
		String[] arg={"/home/hansi/git/mtDNA/mitolib/data/MToolBox/VCF_file.vcf","/home/hansi/git/mtDNA/mitolib/data/MToolBox/VCF_file.txt"};
		VCFReader lr = new VCFReader(arg);
		try {
			lr.build(arg[0], arg[1]);
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
}
