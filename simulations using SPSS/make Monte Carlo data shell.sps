* Encoding: UTF-8.
INPUT PROGRAM.
+ DO REPEAT sim=1 TO 1000.
+    LOOP subj_id=1 TO 50.
+      COMPUTE simstudy_num=sim.
+      COMPUTE group=1.
+      IF  (subj_id > 25) group=2.
+      COMPUTE DV=RV.NORMAL(0,1).
+      END CASE.
+    END LOOP.
+  END REPEAT.
+  END FILE.
END INPUT PROGRAM.
EXECUTE

*an alternative way to do the same thing.
INPUT PROGRAM.
+  LOOP #I = 1 TO 50.
+     DO REPEAT DV = R1 TO R1000.
+          COMPUTE group=1.
+          IF  (#I > 25) group=2.
+          COMPUTE DV=RV.NORMAL(0,1).
+     END REPEAT.                    
+     END CASE.                        
+   END LOOP.
+   END FILE.
END INPUT PROGRAM.
EXECUTE

