"""
File Name:          lr_scheduler_function.py
Project:            bioseq-learning

File Description:

    This file includes all the customized learning rate function used for
    torch.optim.lr_scheduler.LambdaLR

"""
from math import cos, pi


def cosine_annealing_warm_restarts(
        epoch: int,
        eta_max: float = 1.0,
        eta_min: float = 0.1,
        t_0: float = 16.0,
        t_mult: float = 1.002,
        gamma: float = 0.998,
) -> float:
    """lambda function that takes the current number of epoch and returns a
    learning rate based on the cosine annealing with warm restart, proposed
    in https://arxiv.org/abs/1608.03983, with configurable learning rate and
    step size change over the epochs

    reference: https://pytorch.org/docs/stable/optim.html

    :param epoch: current number of epoch, starting from 0
    :type epoch: int
    :param eta_max: maximum learning rate relative to the base learning rate
    :type eta_max: float
    :param eta_min: minimum learning rate relative to the base learning rate
    :type eta_min: float
    :param t_0: initial step size for restart
    :type t_0: float
    :param t_mult: (increase) rate of step size for restart over epochs
    :type t_mult: float
    :param gamma: (decrease) rate of learning rate over epochs
    :type gamma: float
    :return: current learning rate of the epoch
    :rtype: float
    """
    _lr_coef = (gamma ** epoch)
    _t_i = (t_0 * (t_mult ** epoch))
    return eta_min * _lr_coef + 0.5 * _lr_coef * (eta_max - eta_min) * \
        (1 + cos(pi * (epoch % _t_i) / _t_i))
