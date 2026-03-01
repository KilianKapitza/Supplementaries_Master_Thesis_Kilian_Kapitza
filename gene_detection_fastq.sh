#!/bin/bash
# Autor: MitchBioinfo
# Datum: 2025-04-24
# Beschreibung: Verarbeitet mehrere gezippte FastQ-Dateien, filtert sie basierend auf k-mer-Matches
#               und erstellt ein Histogramm der Match-Verteilung über alle Reads.

# ===== KONFIGURATION =====

FASTQ_DIR="./NCF1_20251218_AS_Mat5844_Pat_TP11_fastq_pass_all"        		 			# Ordner mit gezippten FastQ-Dateien
REF_FILE="./NCF1_SNPs_Indels_KK.tsv"        			# Deine Referenzdatei
OUTPUT_FILE="combined_results_NCF1_TH0_KM19_Mat5844_Pat_TP11.fastq" 				# Kombinierte Ausgabedatei
PLOT_FILE="combined_histogram_NCF1_TH0_KM19_Mat5844_Pat_TP11.png"		 			# Ausgabedatei für den Plot
PYTHON_SCRIPT="./gene_detection_multithread_15mer.py" 		# Python-Skript
TEMP_DIR="./temp_results"        					# Temporäres Verzeichnis für Zwischenergebnisse
STATS_DIR="./match_stats"             				# Verzeichnis für Match-Statistiken
THRESHOLD=0                   					# Schwellenwert für die Klassifizierung
PROCESSES=24                    					# Anzahl der zu verwendenden Prozesse


# Erstelle temporäres Verzeichnis, falls es nicht existiert
mkdir -p "$TEMP_DIR"
mkdir -p "$STATS_DIR"

# Leere oder erstelle die kombinierte Ausgabedatei
> "$OUTPUT_FILE"
# Kombinierte Statistikdatei vorbereiten
ALL_STATS="$STATS_DIR/all_matches.stats"
> "$ALL_STATS"

# ===== DATEIEN FINDEN UND ZÄHLEN =====
# Zähle die Anzahl der zu verarbeitenden Dateien
FASTQ_FILES=($(find "$FASTQ_DIR" -name "*.fastq.gz" -o -name "*.fq.gz"))

TOTAL_FILES=${#FASTQ_FILES[@]}

echo "=== FastQ-Dateiverarbeitung ==="
echo "Gefunden: $TOTAL_FILES gezippte FastQ-Dateien im Verzeichnis $FASTQ_DIR"

# Prüfe, ob Dateien gefunden wurden
if [ "$TOTAL_FILES" -eq 0 ]; then
    echo "FEHLER: Keine .fastq.gz oder .fq.gz Dateien im Verzeichnis $FASTQ_DIR gefunden!"
    exit 1
fi

# ===== DATEIEN VERARBEITEN =====
# Initialisiere Zähler
COUNTER=0
PROCESSED_READS=0
FILTERED_READS=0

# Erstelle eine Liste aller temporären Ausgabedateien für die spätere Kombination
TEMP_OUTPUT_FILES=()

echo "=== Starte Verarbeitung ==="
# Verarbeite jede gezippte FastQ-Datei
for FASTQ_FILE in "${FASTQ_FILES[@]}"; do
    COUNTER=$((COUNTER+1))
    FILENAME=$(basename "$FASTQ_FILE")
    TEMP_OUTPUT="$TEMP_DIR/$(basename "$FILENAME" .gz).processed"
    STATS_OUTPUT="$STATS_DIR/$(basename "$FILENAME" .gz).stats"
    TEMP_OUTPUT_FILES+=("$TEMP_OUTPUT")
    
    echo "[$COUNTER/$TOTAL_FILES] Verarbeite $FILENAME..."
    
    # Zähle Anzahl der Reads in der gezippten Datei
    READ_COUNT=$(zcat "$FASTQ_FILE" | awk 'NR%4==1' | wc -l)
    PROCESSED_READS=$((PROCESSED_READS + READ_COUNT))
    echo "  - Datei enthält $READ_COUNT Reads"
    
    # Verarbeite die gezippte Datei mit zcat und sammle Statistiken
    zcat "$FASTQ_FILE" | python "$PYTHON_SCRIPT" --ref "$REF_FILE" --fastq - \
        --output "$TEMP_OUTPUT" --threshold "$THRESHOLD" --processes "$PROCESSES" \
        --stats "$STATS_OUTPUT"
    
    # Prüfe, ob die Verarbeitung erfolgreich war
    if [ $? -ne 0 ]; then
        echo "FEHLER bei der Verarbeitung von $FILENAME! Überspringe diese Datei."
        continue
    fi
    
    # Zähle gefilterte Reads
    if [ -f "$TEMP_OUTPUT" ]; then
        FILTERED_COUNT=$(grep -c "^@" "$TEMP_OUTPUT")
        FILTERED_READS=$((FILTERED_READS + FILTERED_COUNT))
        echo "  - $FILTERED_COUNT Reads nach Filterung behalten"
    fi
    
    # Füge Match-Statistiken zur Gesamtstatistik hinzu
    if [ -f "$STATS_OUTPUT" ]; then
        cat "$STATS_OUTPUT" >> "$ALL_STATS"
    fi
    
    echo "[$COUNTER/$TOTAL_FILES] $FILENAME verarbeitet."
done

echo "=== Verarbeitung abgeschlossen ==="
echo "Insgesamt $PROCESSED_READS Reads verarbeitet, $FILTERED_READS Reads nach Filterung behalten"

# ===== ERGEBNISSE KOMBINIEREN =====
echo "=== Kombiniere gefilterte Reads ==="
# Kombiniere alle temporären Ausgabedateien in eine einzige Datei
for TEMP_FILE in "${TEMP_OUTPUT_FILES[@]}"; do
    if [ -f "$TEMP_FILE" ]; then
        cat "$TEMP_FILE" >> "$OUTPUT_FILE"
    fi
done

echo "Gefilterte Reads wurden in $OUTPUT_FILE kombiniert."

# ===== HISTOGRAMM ERSTELLEN =====
echo "=== Erstelle Histogramm ==="
echo "Erstelle Histogramm aus allen Reads-Statistiken..."

# Prüfe, ob die kombinierte Statistikdatei nicht leer ist
if [ ! -s "$ALL_STATS" ]; then
    echo "FEHLER: Keine Match-Statistiken gefunden! Histogramm kann nicht erstellt werden."
    exit 1
fi

# Erstelle Histogramm mit Python-Inline-Skript
python - <<EOF
import matplotlib.pyplot as plt
import numpy as np

# Lese die Statistiken
with open("$ALL_STATS", 'r') as f:
    matches = [int(line.strip()) for line in f if line.strip()]

threshold = $THRESHOLD
cap_value = 500

# Prüfe, ob Daten vorhanden sind
if not matches:
    print("FEHLER: Keine gültigen Match-Daten gefunden!")
    exit(1)

# Erstelle Histogramm
plt.figure(figsize=(10, 6))
modified_matches = [min(m, cap_value) for m in matches]
exceeding_count = sum(1 for m in matches if m > cap_value)
bins = range(0, cap_value)
counts, _, _ = plt.hist(modified_matches, bins=bins, alpha=0.7, 
                       color='blue', edgecolor='black')

if exceeding_count > 0:
    counts[-2] += exceeding_count
    
plt.yscale('log')
plt.title('Distribution of Reference Matches per Read (All Reads, Log Scale)')
plt.xlabel('Number of Reference Matches')
plt.ylabel('Number of Reads (log scale)')
plt.grid(True, alpha=0.3)
plt.axvline(threshold, color='red', linestyle='dashed', linewidth=2, 
           label=f'Threshold ({threshold})')
plt.legend()
plt.tight_layout()

# Statistiken ausgeben
print(f"Match-Statistiken (alle Reads):")
print(f"  Gesamtzahl Reads: {len(matches)}")
print(f"  Durchschnittliche Matches pro Read: {sum(matches)/len(matches):.2f}")
print(f"  Maximum Matches in einem Read: {max(matches)}")
print(f"  Reads mit >={threshold} Matches: {sum(1 for m in matches if m >= threshold)} ({sum(1 for m in matches if m >= threshold)/len(matches)*100:.2f}%)")
print(f"  Reads mit >{cap_value} Matches: {exceeding_count} ({exceeding_count/len(matches)*100:.2f}%)")

plt.savefig("$PLOT_FILE")
print(f"Histogramm wurde in $PLOT_FILE gespeichert.")
EOF

# ===== AUFRÄUMEN =====
echo "=== Aufräumen ==="
# Frage, ob temporäre Dateien gelöscht werden sollen
read -p "Temporäre Dateien löschen? (j/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Jj]$ ]]; then
    rm -rf "$TEMP_DIR"
    rm -rf "$STATS_DIR"
    echo "Temporäre Dateien wurden gelöscht."
else
    echo "Temporäre Dateien wurden beibehalten."
fi

echo "=== Fertig! ==="
echo "Alle Ergebnisse wurden in $OUTPUT_FILE kombiniert."
echo "Das Histogramm wurde in $PLOT_FILE gespeichert."
