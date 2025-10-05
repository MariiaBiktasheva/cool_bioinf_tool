# Add transcribe
def transcribe(sequence):
    """
    Returns the transcript from the specified sequence
    """
    new_seq = ""

    # Check if this is already a transcript?
    if any(nucl in "Uu" for nucl in sequence):
        print(sequence, "is already a transcript")
        new_seq = sequence

    # Transcribe
    else:
        for nucl in sequence:
            if nucl == "T":
                new_seq = new_seq + "U"
            elif nucl == "t":
                new_seq = new_seq + "u"
            else:
                new_seq = new_seq + nucl
    result = new_seq
    return result
