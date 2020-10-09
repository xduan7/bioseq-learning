"""
File Name:          merge_nni_config.py
Project:            bioseq-learning

File Description:

"""
import logging
from types import MappingProxyType
from typing import Dict, Any


_LOGGER = logging.getLogger(__name__)


def merge_nni_config(
        default_config: MappingProxyType,
        nni_config: Dict[str, Any],
) -> MappingProxyType:
    """merge the default configuration with the NNI configuration returned
    as the next set of hyper-parameters during searching

    :param default_config: default configurations in the format of immutable
    dict (MappingProxyType) of strings to objects; usually generated from a
    *_config.py file
    :type default_config: MappingProxyType
    :param nni_config: next set of hyper-parameters from NNI searching
    :type nni_config: Dict[str, Any]
    :return: immutable dict of configurations of strings to objects
    :rtype: MappingProxyType
    """

    _config_dict: Dict[str, Any] = dict(default_config)

    for _nni_config_key, _nni_config_value in nni_config.items():
        if _nni_config_key not in _config_dict:
            _warning_msg: str = \
                f'NNI configuration with key {_nni_config_key} ' \
                f'is not found in the default configuration. '
            # if the nni config value turned out to be a nested dict
            # recursively call the merge function
            if isinstance(_nni_config_value, dict):
                _warning_msg += f'Trying to add nested dict into config ...'
                _config_dict = dict(merge_nni_config(
                    default_config=MappingProxyType(_config_dict),
                    nni_config=_nni_config_value,
                ))
            else:
                _warning_msg += f'Ignoring ...'
            _LOGGER.warning(_warning_msg)

        else:
            _nni_config_value_type: type = type(_nni_config_value)
            _default_config_value = default_config[_nni_config_key]
            _default_config_value_type: type = type(_default_config_value)

            if _default_config_value_type != _nni_config_value_type:
                _warning_msg: str = \
                    f'NNI configuration with key {_nni_config_key} ' \
                    f'is of type {_nni_config_value}, which is different ' \
                    f'from the type {_default_config_value_type} in the ' \
                    f'default configuration. Performing type casting ...'
                _LOGGER.warning(_warning_msg)

            _config_dict[_nni_config_key] = \
                _default_config_value_type(_nni_config_value)

    return MappingProxyType(_config_dict)
