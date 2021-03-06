"""
File Name:          test_get_module_summary.py
Project:            bioseq-learning

File Description:

"""
import unittest

import torch
import torch.nn as nn
import torch.nn.functional as F

from torchsummary import summary
from utilities.get_module_summary import get_module_summary


_TEST_BATCH_SIZE = 32



class TestGetModuleSummary(unittest.TestCase):
    """unittest class for 'get_module_summary' function
    """
    def test_get_module_summary(self):
        """test 'get_module_summary' function
        """
        class TestMultiInputNet(nn.Module):
            def __init__(self):
                super(TestMultiInputNet, self).__init__()
                self.fc1a = nn.Linear(300, 50)
                self.fc1b = nn.Linear(50, 10)

                self.fc2a = nn.Linear(300, 50)
                self.fc2b = nn.Linear(50, 10)

            def forward(self, x1, x2):
                x1 = F.relu(self.fc1a(x1))
                x1 = self.fc1b(x1)
                x2 = F.relu(self.fc2a(x2))
                x2 = self.fc2b(x2)
                x = torch.cat((x1, x2), -1)
                return F.log_softmax(x, dim=1)

        # model = TestMultiInputNet()
        # input1 = (1, 300)
        # input2 = (1, 300)
        # total_params, trainable_params = summary(
        #     model, [input1, input2], device="cpu")
        # self.assertEqual(total_params, 31120)
        # self.assertEqual(trainable_params, 31120)

        _, summary_str = get_module_summary(
            TestMultiInputNet(),
            (torch.rand(_TEST_BATCH_SIZE, 1, 300).type(torch.float),
             torch.rand(_TEST_BATCH_SIZE, 1, 300).type(torch.float))
        )
        print(summary_str)


if __name__ == '__main__':
    unittest.main()
