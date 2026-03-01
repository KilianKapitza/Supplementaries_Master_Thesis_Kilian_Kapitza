  #!/bin/bash
  
# Pfad des Inputs: bereits durch gene_detection_fastq.sh gefilterte fastq
# Pfad und Benennung des Outputs: Gefilterte Fastq Datei mit erhöht angepasstem Threshold
# Einstellen des Thresholds, also der Anzahl an K-mers die pro Read zu finden sein müssen
# Pfad und Benennung des entstehenden Histogramms
# Aufpassen, dass richtiges conda env aktiv ist! hier ist es: FilterThreshold_Adjustment_env

# Allgemeine Pfade
R_SCRIPT="/home/diablo/Desktop/Kilian_Daten/Testphase/FilterThreshold_Adjustment.r"
OUT_DIR="/home/diablo/Desktop/Kilian_Daten/Testphase/Ergebnisse/NCF1_20251218_AS_Mat5844_Pat_TP11/19-mer/fastq" # Input+Output Path

# Startdatei: TH0
INPUT_FASTQ="${OUT_DIR}/combined_results_NCF1_TH0_KM19_Mat5844_Pat_TP11.fastq"	# Welche Startdatei verwendet werden soll
# Liste der Thresholds, die du testen möchtest
# entweder for THRESHOLD in 50 100 150 200; do für diese 4 Threshold oder for THRESHOLD in $(seq 50 10 200); do für alle Thresholds von 50 bis 200 in 10er Schritten
# #outdated# für schnelleres sorting: erst einmal TH50 und dann $(seq 60 10 200), um Großteil der Reads aussortiert zu haben, wichtig: auch Input dateiname ändern
for THRESHOLD in $(seq 50 10 200); do
  echo ">>> Running threshold ${THRESHOLD}"
  echo "    Input:  ${INPUT_FASTQ}"

  OUTPUT_FASTQ="${OUT_DIR}/combined_results_NCF1_TH${THRESHOLD}_KM19_Mat5844_Pat_TP11.fastq"
  PLOT_FILE="${OUT_DIR}/combined_results_NCF1_TH${THRESHOLD}_KM19_Mat5844_Pat_TP11.png"

  Rscript "$R_SCRIPT" \
    --input "$INPUT_FASTQ" \
    --output "$OUTPUT_FASTQ" \
    --threshold "$THRESHOLD" \
    --plot "$PLOT_FILE"

  echo "    Output: ${OUTPUT_FASTQ}"
  echo ">>> Done: threshold ${THRESHOLD}"
  echo  

# ⚠️ Ab jetzt nächste Runde mit der gerade erzeugten Datei als Input
  INPUT_FASTQ="$OUTPUT_FASTQ"
done


###OUTDATED###
# Nutzung: R_SCRIPT nicht verändern; INPUT_FASTQ pfad verändern je nach Datensatz + wenn TH50 and danach TH60-200 gemacht wird auch INPUT_FASTQ von TH0 zu TH50 ändern
# OUT_DIR in gleichen Ordner wie INPUT_FASTQ; Threshold nach bedarf ändern einzeln mit for THRESHOLD in 50; do oder durchgängig mit for THRESHOLD in $(seq 60 10 200); do
# OUTPUT_FASTQ und PLOT_FILE ebenfalls jedes mal KM und TP ändern


# for THRESHOLD in .... muss nicht mehr geändert werden, da als nächste Startdatei die vorherige Enddatei verwendet wird. Dadurch geht die Geschwindigkeit extrem hinauf
# Namen für OUT_DIR, INPUT_FASTQ, OUTPUT_FASTQ und PLOT_FILE müssen immer noch regelmäßig verändert werden
