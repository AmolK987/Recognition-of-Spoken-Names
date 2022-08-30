#default file 
form Test command line calls
    word File_name amolkr.wav
    real Praat_tstart 0.0
    real Praat_tend   6.0
    real Praat_tstep 0.005
endform


#Read from file: "testamolkr.wav"
#Read from file: "amolkr.wav"
Read from file: file_name$

# Extract the names of the Praat objects
thisSound$ = selected$("Sound")

# Create the Formant Object
select Sound 'thisSound$'
#To Formant (burg)... 0 5 5000 0.025 50
To Formant (burg)... 0 4.5 4500 0.025 50
select Sound 'thisSound$'
# NOTE: To pitch takes only 3 arguments - no advanced arguments allowed
#To Pitch: 0.01, 50, 300
# To Pitch (ac) - allows all the advanced arguments, order is listed below, Help is wrong
pitch_floor = 50
pitch_ceil  = 300
voice_thres = 0.55
#silence_thresh = 0.08
silence_thresh = 0.06
# To Pitch (ac): timestep, pitch floor, max#candidates, very accurate, silence thresh, voice thres, oct cost, oct jmp cst, voice/unv cost, pitch ceil
To Pitch (ac): 0.01, pitch_floor, 15, "off", silence_thresh, voice_thres, 0.01, 0.35, 0.2, pitch_ceil
Rename: "pitch"
select Sound 'thisSound$'
To Intensity: 50, 0.01
Rename: "intensity"

# Create the output file and write the first line.
outputPath$ = "praat_formants.csv"
writeFileLine: "'outputPath$'", "time,F1,F2,F3,F4,pitch,intensity"

#start time
#tstart = 0.09
#tend   = 3.02
#tstep  = 0.00625

#tstart = 0.0
#tend   = 4.9
#tstep  = 0.005

tstart = praat_tstart
tend   = praat_tend
tstep  = praat_tstep

tcurr  = tstart

repeat

   midpoint = tcurr
   tcurr = tcurr + tstep

   # Extract formant measurements
    select Formant 'thisSound$'
    f1 = Get value at time... 1 midpoint Hertz Linear
    f2 = Get value at time... 2 midpoint Hertz Linear
    f3 = Get value at time... 3 midpoint Hertz Linear
    f4 = Get value at time... 4 midpoint Hertz Linear

   # extract pitch
   selectObject: "Pitch pitch"
   pitch = Get value at time: midpoint, "Hertz", "Linear"

   # extract intensity
   selectObject: "Intensity intensity"
   # intensity = Get value at time: midpoint, "Cubic"
   intensity = Get value at time: midpoint, "Linear"

   # Save to a spreadsheet
    appendFileLine: "'outputPath$'", 
                    ...midpoint, ",",
                    ...f1, ",", 
                    ...f2, ",",
                    ...f3, ",", 
                    ...f4, ",",
                    ...fixed$ (pitch, 3), ",",
                    ...fixed$ (intensity, 3)

until (tcurr >= tend)

appendInfoLine: newline$, newline$, "Whoo-hoo! It worked!"