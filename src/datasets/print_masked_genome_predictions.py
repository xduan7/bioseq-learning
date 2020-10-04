"""
File Name:          print_masked_genome_predictions.py
Project:            bioseq-learning

File Description:

"""
import torch
from copy import deepcopy

from colorama import init, Fore

from src.datasets.genome_dataset import INDEX_NUCLEOTIDE_CHAR_DICT


init(autoreset=True)


def print_masked_genome_predictions(
        seq_len: int,
        batch_size: int,
        batch_padded_input: torch.LongTensor,
        batch_padded_target: torch.LongTensor,
        batch_padded_prediction: torch.LongTensor,
):
    _batch_padded_input: torch.LongTensor = \
        deepcopy(batch_padded_input).view(seq_len, batch_size)
    _batch_padded_target: torch.LongTensor = \
        deepcopy(batch_padded_target).view(seq_len, batch_size)
    _batch_padded_prediction: torch.LongTensor = \
        deepcopy(batch_padded_prediction).view(seq_len, batch_size)

    for _i in range(batch_size):

        print('=' * 80)
        # print(f'{_i:3d}-th sample in the batch: ')

        _input = _batch_padded_input[:, _i].view(-1).tolist()
        _target = _batch_padded_target[:, _i].view(-1).tolist()
        _prediction = _batch_padded_prediction[:, _i].view(-1).tolist()

        assert len(_input) == len(_target) == len(_prediction)

        # print predicted sequence with color
        print('pred: ', end='')
        for _j in range(len(_input)):
            _p: str = INDEX_NUCLEOTIDE_CHAR_DICT[_prediction[_j]]
            if _input[_j] == _target[_j]:
                if _prediction[_j] == _target[_j]:
                    print(_p, end='')
                else:
                    print(Fore.YELLOW + _p, end='')
            else:
                if _prediction[_j] == _target[_j]:
                    print(Fore.GREEN + _p, end='')
                else:
                    print(Fore.RED + _p, end='')

        # print targets only if the prediction is wrong
        print('\ntrgt: ', end='')
        for _j in range(len(_input)):
            _t: str = INDEX_NUCLEOTIDE_CHAR_DICT[_target[_j]]
            if _prediction[_j] == _target[_j]:
                print('-', end='')
            else:
                print(_t, end='')
        print('')
