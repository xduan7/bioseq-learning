"""
File Name:          test_merge_nni_config.py
Project:            bioseq-learning

File Description:

"""
import unittest
from typing import Dict, Any
from types import MappingProxyType

from src.utilities import merge_nni_config


class TestMergeNNIConfig(unittest.TestCase):
    """unittest class for 'merge_nni_config' function
    """

    def test_merge_nni_config(self):
        """test 'merge_nni_config' function for torch
        """
        default_config: MappingProxyType = MappingProxyType({
            'config_0': True,
            'config_1': False,
            'config_2': False,
            'config_3': False,
        })
        nni_config: Dict[str, Any] = {
            'ignored_config': -1,
            'config_0': 1,
            'nested_config_1': [
                ['config_1', 1],
                [
                    'nested_config_2',
                    [
                        ['config_2', 1],
                        ['config_3', 1]
                    ]
                ]
            ]
        }
        config = merge_nni_config(
            default_config=default_config,
            nni_config=nni_config,
        )
        print(config)
        assert len(config) == 4
        for _, _v in config.items():
            assert _v


if __name__ == '__main__':
    unittest.main()
