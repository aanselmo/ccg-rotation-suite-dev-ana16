import numpy as np
import excel_write as xl
import pandas as pd
newtry = xl.XLwrite("test_file.xlsx","test_sheet")
s = pd.DataFrame({ 'A' : 1.,
   'B' : pd.Timestamp('20130102'),
   'C' : pd.Series(1,index=list(range(4)),dtype='float32'),
   'D' : np.array([3] * 4,dtype='int32'),
   'E' : pd.Categorical(["test","train","test","train"]),
   'F' : 'foo' })
newtry.writeTable(s.applymap(str))
newtry.newSheet("test_sheet2")
s2 = pd.DataFrame({ 'A' : 1.,
   'B' : pd.Timestamp('20130102'),
   'C' : pd.Series(1,index=list(range(4)),dtype='float32'),
   'D' : np.array([3] * 4,dtype='int32'),
   'E' : pd.Categorical(["test2","train2","test2","train2"]),
   'F' : 'foo' })
newtry.writeTable(s2.applymap(str))
newtry.newSheet("test_sheet3")
s3 = pd.DataFrame({ 'A' : 1.,
   'B' : pd.Timestamp('20130102'),
   'C' : pd.Series(1,index=list(range(4)),dtype='float32'),
   'D' : np.array([3] * 4,dtype='int32'),
   'E' : pd.Categorical(["test3","train3","test3","train3"]),
   'F' : 'foo' })
newtry.writeTable(s3.applymap(str))
newtry.writeTable(s2.applymap(str))
newtry.save()
