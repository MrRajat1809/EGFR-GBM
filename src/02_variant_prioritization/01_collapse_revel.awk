#!/usr/bin/awk -f
BEGIN {
    FS = OFS = "\t"
}

# Skip VEP header comments
/^##/ { next }

# Map header names to column indices for robustness
/^#Uploaded_variation/ {
    for (i = 1; i <= NF; i++) h[$i] = i
    next
}

{
    # 1. Filter for missense variants (includes composite consequences)
    if ($h["Consequence"] !~ /missense_variant/) next

    # 2. Define the unique variant key (Chr, Pos, Alt)
    split($h["Location"], loc, ":")
    chr = loc[1]; pos = loc[2]; alt = $h["Allele"]
    key = chr "\t" pos "\t" alt

    # 3. Extract the numerical REVEL score from the string (e.g., ".,0.257,.,.")
    raw_revel = $h["REVEL_score"]
    revel_val = ""
    split(raw_revel, scores, ",")
    for (i in scores) {
        if (scores[i] ~ /^[0-9.]+$/ && scores[i] != ".") {
            revel_val = scores[i]
            break # Take the first valid numerical score found
        }
    }

    # Skip if no numerical score was found
    if (revel_val == "") next

    # 4. Selection Logic: Prefer Canonical transcript, else take the Highest score
    is_canon = ($h["CANONICAL"] == "YES")

    if (is_canon) {
        best_revel[key] = revel_val
        is_canonical[key] = 1
    } else if (!(key in best_revel)) {
        best_revel[key] = revel_val
    } else if (!is_canonical[key] && revel_val > best_revel[key]) {
        best_revel[key] = revel_val
    }
}

END {
    print "chr", "pos", "alt", "REVEL_score"
    for (k in best_revel) {
        split(k, x, "\t")
        print x[1], x[2], x[3], best_revel[k]
    }
}