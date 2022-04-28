def get_hamming_distance(dna1, dna2):
    return sum([x != y for x, y in zip(dna1, dna2)])


def get_reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[n] for n in dna[::-1]) 


def symbol_to_number(symbol):
    symbols = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    return symbols[symbol]


def number_to_symbol(number):
    return 'ACGT'[number]


def pattern_to_number(pattern):
    if len(pattern) == 0:
        return 0
    return 4 * pattern_to_number(pattern[:-1]) + symbol_to_number(pattern[-1])


def get_quotient(index):
    return index // 4


def get_remainder(index):
    return index % 4


def number_to_pattern(index, k):
    if k == 1:
        return number_to_symbol(index)
    prefix_index = get_quotient(index) 
    symbol = number_to_symbol(get_remainder(index))
    prefix_pattern = number_to_pattern(prefix_index, k-1)
    return prefix_pattern + symbol


def compute_frequencies(dna, k):
    freq_array = [0] * 4**k
    for i in range(len(dna) - k + 1):
        pattern = dna[i:i+k]
        j = pattern_to_number(pattern)
        freq_array[j] += 1
    return freq_array


def get_neighbors(dna, d):
    if d == 0:
        return {dna}
    if len(dna) == 1:
        return {'A','C','G','T'}
    neighborhood = set()
    suffix_neighbors = get_neighbors(dna[1:], d)
    for pattern in suffix_neighbors:
        if get_hamming_distance(dna[1:], pattern) < d:
            for nucleotide in ['A','C','G','T']:
                neighborhood.add(nucleotide+pattern)
        else:
            neighborhood.add(dna[1]+pattern)
    return neighborhood


def get_approximate_pattern_count(seq, pattern, d):
    count = 0
    for i in range(len(seq) - len(pattern) + 1):
        if get_hamming_distance(pattern, seq[i: i + len(pattern)]) <= d:
            count += 1 
    return count


