package genepi.mitolib.haplogrep;

import java.io.IOException;
import java.util.HashMap;

import genepi.base.Tool;
import genepi.io.table.reader.CsvTableReader;
import genepi.io.text.LineWriter;

public class VariantsToHsd extends Tool {

	public VariantsToHsd(String[] args) {
		super(args);
	}

	@Override
	public void init() {
		System.out.println("Create Haplogrep Input for stable variants");
	}

	@Override
	public void createParameters() {

		addParameter("in", "variants.txt");
		addParameter("out", "output file for HaploGrep");
	}

	@Override
	public int run() {

		String in = (String) getValue("in");

		String out = (String) getValue("out");

		CsvTableReader reader = new CsvTableReader(in, '\t');

		HashMap<String, String> profiles = new HashMap<String, String>();

		while (reader.next()) {

			double level = reader.getDouble("Variant-Level");

			// only take stable positions!
			if (level > 0.5) {

				String id = reader.getString("SampleID");

				int pos = reader.getInteger("Pos");

				String variant = reader.getString("Variant");

				String snp = pos + variant;

				if (profiles.get(id) == null) {

					profiles.put(id, snp);
					
				} else {

					String value = profiles.get(id);

					profiles.put(id, value + "\t" + snp);

				}

			}
		}

		LineWriter writer = null;
		try {
			writer = new LineWriter(out);

			for (String key : profiles.keySet()) {

				StringBuilder hsdLine = new StringBuilder();

				hsdLine.append(key + "\t" + "1-16569" + "\t" + "?" + "\t");

				hsdLine.append(profiles.get(key));

				writer.write(hsdLine.toString());
			}

			writer.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.out.println("Done. Written to "+out);
		
		return 0;

	}

	public static void main(String[] args) throws IOException {

		// new AroFilter2(args).start();
		VariantsToHsd test = new VariantsToHsd(
				new String[] { "--in", "data/haplogrep/variants.txt", "--out", "data/haplogrep/samples.hsd" });
		test.start();

	}

}
