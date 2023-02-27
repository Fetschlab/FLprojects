%{
SJ 04/2022
updated 06/2022
Notes on the neural data processing pipeline for the Ripple Grapevine/MP
rig, explaining the logic and workflow.

%============
1. RECORDING AND RAW DATA

Neural data and timestamps of key behavioral events are recorded on Trellis
system. 

Timestamps and online sorted digital spike events are stored in .nev files
1kHz analog data (eye signals and IMU accelerometer signals) are in .ns2
30kHz neural data (raw broadband, lfp, and high-pass filtered) are in .ns5

These raw data files should be moved to the date folder created by Trellis, and then copied to NAS, immediately after recording
is complete. Note that time starts at zero at the start of each Trellis recording.

%============
2. PRE-PROCESSING

1) input recording notes into a cell block within createTrellisInfo.m for
the given recording day (see template), and run the cell. This will create and save an info
struct containing all the header info/metadata for the day's session.

2) processTrellisData will then be run using this info file, to create an
nsEvents struct (saved in [subject][date]dots3DMP[filenum]_RippleEvents.mat), and a binary int16 file for
Kilosort (saved as [subject][date]_[recording set number] (if createBinaryFiles == 1)
This will create one binary file per recording 'set'/block (from one or more Ripple Events files), so that
waveforms from the same recording depth across different PLDAPS experiments
can be sorted together, saved as _1, _2 etc. See *SPIKE SORTING* note
below. It is therefore important that rec_group is correctly specified in
info file so that pipeline knows which files belong together.

Notes on nsEvents (saved in *_RippleEvents.mat file)
- contains trial condition and outcome information, and timing of key
behavioral events (all in seconds in Trellis recording time)

since trial conditions are sent as binary codes, all except modality (i.e. heading, delta, coherence) are
coded as the index in condition list, rather than the actual value. Offline recovery of actual condition
value is done using nsEventConditions.m during post-processing (taken care
of in 3b.)

- for online window discriminated spikes (if any), nsEvents will contain a spkData field with
spike timing and cluster information (if chanInterest field in info is not empty)

- for offline alignment purposes (if multiple Trellis recordings were made at the same depth i.e. in the same set/group), nsEvents will contain a field called
analogData, with important header information about the neural data. For
aligning spikes to behavioral events, spike times extracted from Kilosort need to be corrected by the
'startTime' in the relevant nsEvents.analogData for aligning to events in
nsEvents, or the nsEvents can be corrected by the timestamps in order to
show data from multiple recordings together. See *SPIKE SORTING* note below

%============
3. NEXT STEP PROCESSING

a. Offline spike sorting for probe recordings should be done with Kilosort,
pointing to the relevant binary file for a given date and recording set.
(main_kilosort.m, SJ parameters on local Trellis computer, with appropriate config files set up)
After manual curation, necessary .phy output files will be in the _[set]
folder where the binary file is located (on NAS).

For single electrode recordings, .ns5 files can be loaded into MKsort in
Matlab, and sorted manually.

b. creating neural data structure, on any computer
for a given subject, dateRange, and paradigm, run
dots3DMP_NeuralPreProcessing, which calls getNeuralEventsInfo and createSessionData to
generate a cell containing neural info. See these codes for more info on
neural data organization. 
Briefly, getNeuralEventsInfo locally downloads the info and RippleEvents
files, createSessionData uses these and the sorted spike files (e.g .npy from
Kilosort, waveforms_[ch].mat from MKsort) to create a struct containing neural Data and metadata for each
recording session/set.

*** NOTE ON SPIKE SORTING ***

Generally it is better for Kilosort to handle recordings done at the same
depth (with presumably the same neurons) as one file, instead of a user
trying to manually reconcile clusters across separately sorted files.
a) it will require a lot of effort and be error-prone to try and reconcile
cluster ids across separately sorted recordings (if the same units are
supposed to be in both recordings). b) longer recordings provide kilosort
with more data to work with for templates etc.

On the other hand, it can be more straightforward to make sure each Trellis
recording is linked to just one PLDAPS file on the experimental computer.
(This means stopping and starting Trellis recordings for each paradigm e.g.
tuning and task, even at the same depth).
(Why is it more straightforward?)
- 1. Because we need to reconcile the PLDAPS conditions list with the
nsEvents one (and also nexonar potentially)
- 2. The PLDAPS events are different for different paradigms, so
constructing nsEvents is more parsimonious / cleaner
- 3. We can decide offline how to combine recordings for
analysis/visualization e.g. two dots3DMP recordings for different modalities

The alternative is to just run the Trellis recording for the whole duration of recording different 'paradigms' at the same location, and deal with the
events offline. In this case, pipeline will still create one binary file for a given recording set for
Kilosort purposes, and for creating the neural data cell we just sub-select which segments of task events and spike times we are interested in.

As of May 2022, the typical procedure is to run one long Trellis
recording. nsEvents creation can now handle multiple paradigms. Binary file creation for Kilosort requires sequential writing to
file of the .ns5 raw data because the files are too large to load into
Matlab memory. Post-processing in createSessionData sub-selects the task
events and roughly associated spike times for one paradigm to create the
fields for different paradigms in neural dataCell.

%}
