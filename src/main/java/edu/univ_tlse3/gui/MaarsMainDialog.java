package edu.univ_tlse3.gui;

import edu.univ_tlse3.cellstateanalysis.SetOfCells;
import edu.univ_tlse3.display.SOCVisualizer;
import edu.univ_tlse3.maars.MAARS;
import edu.univ_tlse3.maars.MAARSNoAcq;
import edu.univ_tlse3.maars.MaarsParameters;
import edu.univ_tlse3.utils.FileUtils;
import edu.univ_tlse3.utils.GuiUtils;
import edu.univ_tlse3.utils.IOUtils;
import ij.IJ;
import ij.gui.YesNoCancelDialog;
import ij.plugin.frame.RoiManager;
import mmcorej.CMMCore;
import org.apache.commons.math3.util.FastMath;
import org.micromanager.internal.MMStudio;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.IOException;
import java.util.Map;
import java.util.concurrent.*;

/**
 * Class to create and display a dialog to get parameters of the whole analysis
 *
 * @author Tong LI
 */

public class MaarsMainDialog extends JFrame implements ActionListener {

    private final MMStudio mm;
    private final CMMCore mmc;
    private MaarsParameters parameters;
    private JButton okMainDialogButton;
    private JButton showDataVisualizer_;
    private JButton segmButton;
    private JButton fluoAnalysisButton;
    private JFormattedTextField savePathTf;
    private JFormattedTextField fluoAcqDurationTf;
    private JCheckBox postAnalysisChk_;
    private JRadioButton dynamicOpt;
    private JRadioButton staticOpt;
    private MaarsFluoAnalysisDialog fluoDialog_;
    private MaarsSegmentationDialog segDialog_;
    private SetOfCells soc_ = new SetOfCells();
    private SOCVisualizer socVisualizer_;
    private JButton stopButton_;
    private MAARS maars_;
    private MAARSNoAcq maarsNoAcq_;
    private CopyOnWriteArrayList<Map<String, Future>> tasksSet_ = new CopyOnWriteArrayList<>();
    private JCheckBox saveParametersChk_;

    /**
     * Constructor
     *
     * @param mm         : graphical user interface of Micro-Manager
     * @param parameters :MaarsParameters
     */
    public MaarsMainDialog(MMStudio mm, MaarsParameters parameters) {
        super("Mitosis Analysing And Recording System - MAARS");
        // ------------initialization of parameters---------------//

        this.mm = mm;
        this.mmc = mm.core();
        this.parameters = parameters;

        IJ.log("create main dialog ...");
        setDefaultLookAndFeelDecorated(true);
        setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

        // set minimal dimension of mainDialog

        int maxDialogWidth = 350;
        int maxDialogHeight = 600;
        Dimension minimumSize = new Dimension(maxDialogWidth, maxDialogHeight);
        setMinimumSize(minimumSize);

        // Exploration Label

        JPanel multiPositionPanel = new JPanel(new GridLayout(2, 1));
        multiPositionPanel.setBackground(GuiUtils.bgColor);
        multiPositionPanel.setBorder(GuiUtils.addPanelTitle("Path to position list (.pos)"));

        JFormattedTextField posListTf = new JFormattedTextField(String.class);
        posListTf.setText(parameters.getPathToPositionList());
        multiPositionPanel.add(posListTf);

        JPanel posListActionPanel = new JPanel(new GridLayout(1, 0));
        final JButton editPositionListButton = new JButton("Generate...");
        editPositionListButton.addActionListener(e -> mm.showPositionList());
        posListActionPanel.add(editPositionListButton);

        multiPositionPanel.add(posListActionPanel);
        // analysis parameters label

        JPanel analysisParamPanel = new JPanel(new GridLayout(2, 1));
        analysisParamPanel.setBackground(GuiUtils.bgColor);
        analysisParamPanel.setBorder(GuiUtils.addPanelTitle("Analysis parameters"));
        analysisParamPanel.setToolTipText("Set parameters");

        // segmentation button

        JPanel segPanel = new JPanel(new GridLayout(1, 0));
        segmButton = new JButton("Segmentation");
        segmButton.addActionListener(this);
        segPanel.add(segmButton);
        analysisParamPanel.add(segPanel);

        // fluo analysis button

        JPanel fluoAnalysisPanel = new JPanel(new GridLayout(1, 0));
        fluoAnalysisButton = new JButton("Fluorescence analysis");
        fluoAnalysisButton.addActionListener(this);
        fluoAnalysisPanel.add(fluoAnalysisButton);
        analysisParamPanel.add(fluoAnalysisPanel);

        // strategy panel (2 radio button + 1 textfield + 1 label)

        JPanel strategyPanel = new JPanel(new GridLayout(1, 0));
        strategyPanel.setBorder(GuiUtils.addPanelTitle("Strategy"));
        strategyPanel.setToolTipText("Which strategy to use");

//      strategyPanel.setBackground(panelColor);
        dynamicOpt = new JRadioButton("Dynamic");
        dynamicOpt.setSelected(parameters.useDynamic());
        staticOpt = new JRadioButton("Static");
        staticOpt.setSelected(!parameters.useDynamic());

        dynamicOpt.addActionListener(this);
        staticOpt.addActionListener(this);

        ButtonGroup group = new ButtonGroup();
        group.add(dynamicOpt);
        group.add(staticOpt);

        strategyPanel.add(staticOpt);
        strategyPanel.add(dynamicOpt);
        fluoAcqDurationTf = new JFormattedTextField(Double.class);
        fluoAcqDurationTf.setValue(parameters.getFluoParameter(MaarsParameters.TIME_LIMIT));
        strategyPanel.add(fluoAcqDurationTf);
        strategyPanel.add(new JLabel("min", SwingConstants.CENTER));
        strategyPanel.setBackground(GuiUtils.bgColor);

        // checkbox : update or not MAARS parameters

        JPanel chkPanel = new JPanel(new GridLayout(1, 0));
        chkPanel.setBackground(GuiUtils.bgColor);
        chkPanel.setBorder(GuiUtils.addPanelTitle("Options"));
        chkPanel.setToolTipText("check post analysis if without microscope, check parameter if don't want to replace last parameters");

        saveParametersChk_ = new JCheckBox("Save parameters", true);
        postAnalysisChk_ = new JCheckBox("Post analysis", false);
        postAnalysisChk_.addActionListener(this);
        chkPanel.add(postAnalysisChk_);
        chkPanel.add(saveParametersChk_);

        // Saving path Panel

        JPanel savePathPanel = new JPanel(new GridLayout(1, 0));
        savePathPanel.setBackground(GuiUtils.bgColor);
        savePathPanel.setBorder(GuiUtils.addPanelTitle("Acquisition root folder :"));
        savePathPanel.setToolTipText("Path to saving folder");

        // Saving Path textfield

        savePathTf = new JFormattedTextField(parameters.getSavingPath());
        savePathTf.setMaximumSize(new Dimension(maxDialogWidth, 1));
        savePathPanel.add(savePathTf);

        // show visualiwer acquisitions button

        JPanel stopAndVisualizerButtonPanel_ = new JPanel(new GridLayout(1, 0));
        stopAndVisualizerButtonPanel_.setBackground(GuiUtils.bgColor);
        stopAndVisualizerButtonPanel_.setBorder(GuiUtils.addPanelTitle("Visualizer and stop"));
        showDataVisualizer_ = new JButton("Show visualizer");
        showDataVisualizer_.addActionListener(this);
        stopAndVisualizerButtonPanel_.add(showDataVisualizer_);

        //

        stopButton_ = new JButton("Stop");
        stopButton_.addActionListener(this);
        stopAndVisualizerButtonPanel_.add(stopButton_);


        // Ok button to run

        JPanel okPanel = new JPanel(new GridLayout(1, 0));
        okMainDialogButton = new JButton("Go !");
        okMainDialogButton.setBackground(GuiUtils.butColor);
        okMainDialogButton.addActionListener(this);
        getRootPane().setDefaultButton(okMainDialogButton);
        okPanel.add(okMainDialogButton);

        // ------------set up and add components to Panel then to Frame---------------//

        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.add(multiPositionPanel);
        mainPanel.add(analysisParamPanel);
        mainPanel.add(strategyPanel);
        mainPanel.add(chkPanel);
        mainPanel.add(savePathPanel);
        mainPanel.add(okPanel);
        mainPanel.add(stopAndVisualizerButtonPanel_);
        add(mainPanel);
        IJ.log("Done.");
        pack();
    }

    public static void waitAllTaskToFinish(CopyOnWriteArrayList<Map<String, Future>> tasksSet) {
        for (Map<String, Future> aFutureSet : tasksSet) {
            for (String channel : aFutureSet.keySet()) {
                try {
                    aFutureSet.get(channel).get();
                } catch (InterruptedException | ExecutionException e) {
                    IOUtils.printErrorToIJLog(e);
                }
                IJ.showStatus("Terminating analysis...");
            }
        }
        IJ.log("Spot detection finished! Proceed to saving and analysis...");
    }

    /**
     * @return graphical user interface of Micro-Manager
     */
    private MMStudio getMM() {
        return mm;
    }

    /**
     * method to save the parameters entered
     */
    private void saveParameters() {
        if (!savePathTf.getText().equals(parameters.getSavingPath())) {
            parameters.setSavingPath(savePathTf.getText());
        }
        if (!fluoAcqDurationTf.getText().equals(parameters.getFluoParameter(MaarsParameters.TIME_LIMIT))) {
            parameters.setFluoParameter(MaarsParameters.TIME_LIMIT, fluoAcqDurationTf.getText());
        }
        try {
            if (saveParametersChk_.isSelected()) {
                parameters.save();
            }
            if (FileUtils.exists(parameters.getSavingPath())) {
                parameters.save(parameters.getSavingPath());
            }
        } catch (IOException e) {
            IJ.error("Could not save parameters");
        }
    }

    /**
     * method to set the strategy selected
     */
    private void setAnalysisStrategy() {

        if (dynamicOpt.isSelected()) {
            parameters.setFluoParameter(MaarsParameters.DYNAMIC, "" + true);
        } else if (staticOpt.isSelected()) {
            parameters.setFluoParameter(MaarsParameters.DYNAMIC, "" + false);
        }
    }

    private int overWrite(String path) {
        int overWrite = 0;
        if (FileUtils.exists(path + File.separator + "X0_Y0" + File.separator + "MMStack.ome.tif")) {
            overWrite = JOptionPane.showConfirmDialog(this, "Overwrite existing acquisitions?");
        }
        return overWrite;
    }

    private SOCVisualizer createVisualizer() {
        final SOCVisualizer socVisualizer = new SOCVisualizer();
        socVisualizer.createGUI(soc_);
        return socVisualizer;
    }

    private void setSkipTheRest(Boolean stop) {
        if (maarsNoAcq_ != null) {
            maarsNoAcq_.stop_.set(stop);
        }
        if (maars_ != null) {
            maars_.skipAllRestFrames = stop;
        }
    }

    @Override
    public void actionPerformed(ActionEvent e) {
        if (e.getSource() == okMainDialogButton) {
            if (socVisualizer_ == null) {
                socVisualizer_ = createVisualizer();
                if (parameters.useDynamic()) {
                    socVisualizer_.setVisible(true);
                }
            }
            saveParameters();
            setSkipTheRest(false);
            if (postAnalysisChk_.isSelected()) {
                ExecutorService es = Executors.newSingleThreadExecutor();
                es.execute(new MAARSNoAcq(parameters, socVisualizer_, soc_));
                es.shutdown();
            } else {
                if (overWrite(parameters.getSavingPath()) == JOptionPane.YES_OPTION) {
                    ExecutorService es = Executors.newSingleThreadExecutor();
                    es.execute(new MAARS(mm, mmc, parameters, socVisualizer_, tasksSet_, soc_));
                    es.shutdown();
                }
            }
        } else if (e.getSource() == segmButton) {
            saveParameters();
            if (segDialog_ != null) {
                segDialog_.setVisible(true);
            } else {
                segDialog_ = new MaarsSegmentationDialog(this, parameters, mm);
            }

        } else if (e.getSource() == fluoAnalysisButton) {
            saveParameters();
            if (fluoDialog_ != null) {
                fluoDialog_.setVisible(true);
            } else {
                fluoDialog_ = new MaarsFluoAnalysisDialog(this, mm, parameters);
            }
        } else if (e.getSource() == dynamicOpt) {
            setAnalysisStrategy();
            fluoAcqDurationTf.setEditable(true);
        } else if (e.getSource() == staticOpt) {
            setAnalysisStrategy();
            fluoAcqDurationTf.setEditable(false);
        } else if (e.getSource() == showDataVisualizer_) {
            if (socVisualizer_ == null) {
                socVisualizer_ = createVisualizer();
            }
            socVisualizer_.setVisible(true);
        } else if (e.getSource() == stopButton_) {
            YesNoCancelDialog yesNoCancelDialog = new YesNoCancelDialog(this, "Abandon current acquisition?",
                    "Stop current analysis ?");
            yesNoCancelDialog.setAlwaysOnTop(true);
            if (yesNoCancelDialog.yesPressed()) {
                setSkipTheRest(true);
                RoiManager roiManager = RoiManager.getInstance();
                roiManager.runCommand("Select All");
                roiManager.runCommand("Delete");
                roiManager.reset();
                roiManager.close();
                soc_.reset();
                socVisualizer_.cleanUp();
                socVisualizer_.setVisible(false);
                socVisualizer_.createGUI(soc_);
            }
        } else {
            IJ.log("MAARS don't understand what you want, sorry");
        }
    }
}
