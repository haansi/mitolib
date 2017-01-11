package gui;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;

import javax.swing.BoxLayout;
import javax.swing.GroupLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.SpringLayout;
import javax.swing.SwingUtilities;

import genepi.bam.VariantBuilder;

public class greenGUI extends JFrame {

	public greenGUI() {
		super("greenVC-GUI");

		setLayout(new GridLayout(0, 2));

		// set up a file picker component
		JLabel bamLabel = new JLabel("Please choose a bam/sam file");
		final JFilePicker bamFile = new JFilePicker("Browse...");
		bamFile.setMode(JFilePicker.MODE_OPEN);
		bamFile.addFileTypeFilter(".bam", "BAM File");
		bamFile.addFileTypeFilter(".sam", "SAM File");
		bamLabel.setHorizontalAlignment(JLabel.CENTER);
		add(bamLabel);
		add(bamFile);

		JLabel fastaLabel = new JLabel("Please choose a fasta file");
		final JFilePicker fastaFile = new JFilePicker("Browse...");
		fastaFile.setMode(JFilePicker.MODE_OPEN);
		fastaFile.addFileTypeFilter(".fasta", "FASTA File");
		fastaLabel.setHorizontalAlignment(JLabel.CENTER);
		add(fastaLabel);
		add(fastaFile);

		JLabel outputLabel = new JLabel("Please choose the output folder");
		final JFilePicker outPutFolder = new JFilePicker("Browse...");
		outPutFolder.setMode(JFilePicker.MODE_SAVE_FOLDER);
		outputLabel.setHorizontalAlignment(JLabel.CENTER);
		add(outputLabel);
		add(outPutFolder);

		JPanel help = new JPanel(new SpringLayout());
		JLabel labelVAF = new JLabel("VAF (Variant Allele Frequency - heteroplasmy level):");
		final JTextField textFieldVAF = new JTextField(5);
		textFieldVAF.setText("0.05");
		labelVAF.setHorizontalAlignment(JLabel.CENTER);
		add(labelVAF);
		help.add(textFieldVAF);
		add(help);

		JPanel help2 = new JPanel(new SpringLayout());
		JLabel labelQual = new JLabel("Quality (Phred Score):");
		final JTextField textFielQual = new JTextField(2);
		textFielQual.setText("20");
		labelQual.setHorizontalAlignment(JLabel.CENTER);
		help2.add(textFielQual);
		add(labelQual);
		add(help2);

		JLabel runlabel = new JLabel("");
		JButton runButton = new JButton("Run");

		runButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				if (textFielQual.getText().length() > 0 && textFieldVAF.getText().length() > 0
						&& !bamFile.getSelectedFilePath().isEmpty() && !fastaFile.getSelectedFilePath().isEmpty() && !outPutFolder.getSelectedFilePath().isEmpty()) {
					String in = (String) bamFile.getSelectedFilePath();
					String out = (String) outPutFolder.getSelectedFilePath();
					String ref = (String) fastaFile.getSelectedFilePath();
				
					double vaf = Double.valueOf(textFieldVAF.getText()) ;
					int qual = Integer.valueOf(textFielQual.getText());

					VariantBuilder builder = new VariantBuilder(in);
					builder.setReference(ref);
					builder.setOutDirectory(out);
					builder.setVaf(vaf);
					try {
						 builder.build();
						 JOptionPane.showMessageDialog(new JFrame(), "Data processed", "Info",
									JOptionPane.INFORMATION_MESSAGE);
					} catch (MalformedURLException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
						 JOptionPane.showMessageDialog(new JFrame(), "Data not processed", "Error",
									JOptionPane.ERROR_MESSAGE);
					} catch (IOException e1) {
						// TODO Auto-generated catch block
						JOptionPane.showMessageDialog(new JFrame(), "Data not processed", "Error",
								JOptionPane.ERROR_MESSAGE);
						e1.printStackTrace();
					
					}
					
				} else {
					JOptionPane.showMessageDialog(new JFrame(), "Not all data provided", "Attention",
							JOptionPane.WARNING_MESSAGE);

				}
			}

		});

		add(runlabel);
		add(runButton);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(920, 200);
		setLocationRelativeTo(null); // center on screen
	}

	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				new greenGUI().setVisible(true);
			}
		});
	}

}