# TrillOMatic
This is a matlab implementation of the trillOmatic: 
a software for automatic segmentation and analysis of bird song trills

## Running the trillOmatic
first, run the matlab script `addpath_birds.m`. This will automatically add all relevant files 
to the system path. Afterwards, type `trillOmatic` to start.

## Working with the trillOmatic

### File Selection
The most straightforward way to load a file in trillOmatic is by pressing the `Other...` push button,
and selecting the relevant file from the dialog box that opens.
You can also add specific files for quick run to the quick selection box. To do this, first, make 
sure that the files you wish to analyze are under the `trillOmatic_path/sounds/` sub folder.
Then, open the file `trillOmatic_path/birdsoft/filenames.csv` and add the file name in the last column.

### Main Windows

The middle window contains the waveform of the audio file you selected. The bottom window contains
its spectrogram, and the top window shows its Teager energy. 

### Segmentation Algorithms

There are two automatic segmentation algorithms the trillOmatic currently offers. They are fully 
described in Hagai Barmatz's Msc thesis, which can be found here: 
`https://www.tau.ac.il/~stoledo/Theses/Hagai_Master_Thesis_4_0.pdf`

### Manual Segmenatation
In some cases, the user may want to refine the automatic segmentation produced by the algorithm.
This can be done by simply interacting with the waveform window.

#### Adding a New Syllable
By simply clicking on the waveform window, where no segmentation is present, a new pair of 
segmentation indicators will be introduced.
#### Deleting an Existing Syllable
By clicking on one of the indicators of an existing syllable, a dialog will be presented, asking 
to confirm deleting the segmentation.
#### Editing a Segmentation
By clicking on a segmentation indicator and dragging it, the segmentation can be manually adjusted.

### Estimating parameters
By clicking on the `Export Data` push button, a table with parameters of interest will be presented.
Exporting the information to an `.xlsx` file for further analysis is also possible.

For each segmented syllable, the following parameters are presented:
Start time,	End Time, Max f0,	Min f0,	Bandwidth,	f slope,	Duration,	Interval Duration,	Period Time

In addition, the following parameters are calculated for the whole trill:
Syllable Count,	Trill Duration,	Syllable Rate,	Syllable Duration,	Interval Duration,	PRI,	Syllable Bandwidth,	Frequency Slope (Syllable),	Max Amplitude,	Attack Time,	Attack Relative Time,	Release Time,	Release Relative Time,	Max Band,	Min Band,	f0 Start Max,	f0 End Max,	Delta Max,	Max Slope,	f0 Start Min,	f0 End Min,	Delta Min,	Min Slope

### Band Pass Filter
In some cases, irrelevant birds and other noise may be present in the recording. If this interferes
with coherent automatic segmentation, noise can be filtered out using a band pass filter. 
Set the appropriate high and low frequencies (in kHz) and check the `BP Filter` checkbox. After 
the signal is filtered to your satisfaction, run the segmentation algorithms and obtain the new 
segmentation.

### Spectrogram options
Spectrogram paramters such as window type and size can be adjusted in the dialog that opens by 
clicking the `Spect. Options` push button.

### Baseline
By default, the baseline is filtered out from the waveform. If it contains important information,
one can show the waveform with the original baseline by checking the `Baseline` checkbox.

### Sound
By pressing the `Sound` push button, the audio file will be played through your system's speakers.