import argparse
from collections import Counter

# TODO: Adapt Theseus to this format?
def read_msa(filename):
  """Read MSA file (FASTA format)."""
  sequences = {}
  current_id = None
  with open(filename, 'r') as f:
    for line in f:
      line = line.strip()
      if line.startswith('>'):
        current_id = line[1:]
        sequences[current_id] = ""
      elif current_id is not None:
        sequences[current_id] += line
  return sequences


# Define scoring used by Sum of Pairs
sp_m      =  1
sp_x      = -1
sp_g      = -1
sp_g_to_g =  0

def sum_of_pairs(sequences):
  """Calculate Sum of Pairs score."""
  seq_list = list(sequences.values())
  if not seq_list or len(seq_list[0]) == 0:
    return 0

  score = 0
  align_len = len(seq_list[0])

  # Sum of pairs score is computed over all pairs of nucleotides or gaps per column
  for col in range(align_len):
    for i in range(len(seq_list)):
      for j in range(i + 1, len(seq_list)):
        # Match or mismatch
        if seq_list[i][col] != '-' and seq_list[j][col] != '-':
          if seq_list[i][col] == seq_list[j][col]:
            score += sp_m
          else:
            score += sp_x
        # Gap or gap-to-gap
        else:
          if seq_list[i][col] != seq_list[j][col]:
            score += sp_g
          else:
            score += sp_g_to_g

  return score


def scaled_sum_of_pairs(sequences):
  # TODO: Which length do we use to normalize
  first_seq = list(sequences.values())[0]
  return sum_of_pairs(sequences)/len(first_seq)


# Counts the total number of fully conserved columns. That is, how many columns have the same
# non-gap character for all sequences.
def total_column(sequences):
  """Calculate Total Column score (fully conserved columns)."""
  seq_list = list(sequences.values())
  if not seq_list or len(seq_list[0]) == 0:
    return 0

  score = 0
  align_len = len(seq_list[0])
  for col in range(align_len):
    column = [seq[col] for seq in seq_list]
    if len(set(column)) == 1 and column[0] != '-':
      score += 1
  return score

# Compare to a reference sequence
# def q_score(sequences):
#   """Calculate Q-score (fraction of identical residues in columns)."""
#   seq_list = list(sequences.values())
#   if not seq_list or len(seq_list[0]) == 0:
#     return 0

#   n = len(seq_list)
#   align_len = len(seq_list[0])
#   total_pairs = n * (n - 1) // 2 * align_len

#   if total_pairs == 0:
#     return 0

#   identical_pairs = sum_of_pairs(sequences)
#   return identical_pairs / total_pairs if total_pairs > 0 else 0


def main():
  parser = argparse.ArgumentParser(description="Compute MSA accuracy metrics")
  parser.add_argument("msa_file", help="Path to MSA file (FASTA format)")
  args = parser.parse_args()

  # Read sequences
  sequences = read_msa(args.msa_file)

  # Compute accuracy values
  scaled_sp = scaled_sum_of_pairs(sequences)
  sp        = sum_of_pairs(sequences)
  tc        = total_column(sequences)
  # q = q_score(sequences)

  # Print values
  print(f"Scaled Sum of Pairs (SP): {scaled_sp:.4f}")
  print(f"Sum of Pairs (SP): {sp}")
  print(f"Total Column (TC): {tc}")
  # print(f"Q-score: {q:.4f}")

if __name__ == "__main__":
  main()