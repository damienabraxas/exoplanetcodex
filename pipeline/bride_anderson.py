"""Experimental Anderson acceleration of the physical direct-iteration map."""
import numpy as np


class LogPopulationAnderson:
    """Accelerate log fractions; keep positivity and particle conservation explicit.

    History contains actual direct-map input/output pairs, never approximate MALI
    equations. Rejected proposals return the unaltered direct-map output. Final
    physical-rate/radiation qualification is the caller's responsibility.
    """
    def __init__(self, *, memory, delay, rcond, max_log_step, step_policy='reject'):
        if memory<1 or delay<0 or not 0<rcond<1 or max_log_step<=0:
            raise ValueError('Invalid Anderson configuration')
        self.memory=memory;self.delay=delay;self.rcond=rcond
        self.max_log_step=max_log_step;self.history=[];self.calls=0
        if step_policy not in ['reject','backtrack']:raise ValueError('Unknown Anderson step policy')
        self.step_policy=step_policy

    def update(self, previous, direct, total):
        previous=np.asarray(previous);direct=np.asarray(direct);total=np.asarray(total)
        if previous.shape!=direct.shape or previous.ndim!=2 or total.shape!=(direct.shape[1],):
            raise ValueError('Inconsistent Anderson state')
        if any(not np.isfinite(v).all()or np.any(v<=0)for v in [previous,direct,total]):
            raise ValueError('Invalid Anderson population')
        x=np.log(previous/total);g=np.log(direct/total);f=g-x
        self.calls+=1;self.history.append((g.copy(),f.copy()))
        self.history=self.history[-(self.memory+1):]
        info=dict(attempted=False,accepted=False,rank=None,proposal_max_log_step=None)
        if self.calls<=self.delay or len(self.history)<2:return direct.copy(),info
        info['attempted']=True
        dg=np.column_stack([(b[0]-a[0]).ravel()for a,b in zip(self.history[:-1],self.history[1:])])
        df=np.column_stack([(b[1]-a[1]).ravel()for a,b in zip(self.history[:-1],self.history[1:])])
        coefficients,_,rank,_=np.linalg.lstsq(df,f.ravel(),rcond=self.rcond)
        info['rank']=int(rank)
        proposal=g-(dg@coefficients).reshape(g.shape)
        step=float(np.max(abs(proposal-x)));info['proposal_max_log_step']=step if np.isfinite(step)else None
        if not np.isfinite(proposal).all():
            self.history=self.history[-1:];return direct.copy(),info
        if self.step_policy=='backtrack':
            direction=proposal-x;factor=1.
            for _ in range(61):
                trial=x+factor*direction;trial-=np.logaddexp.reduce(trial,axis=0)
                actual_step=float(np.max(abs(trial-x)))
                if np.isfinite(trial).all()and actual_step<=self.max_log_step:
                    proposal=trial;info['backtracking_factor']=factor;info['normalized_max_log_step']=actual_step;break
                factor*=.5
            else:
                self.history=self.history[-1:];return direct.copy(),info
        else:
            if step>self.max_log_step:
                self.history=self.history[-1:];return direct.copy(),info
            proposal-=np.logaddexp.reduce(proposal,axis=0)
        candidate=np.exp(proposal)*total
        if not np.isfinite(candidate).all()or np.any(candidate<=0):
            self.history=self.history[-1:];return direct.copy(),info
        info['accepted']=True
        return candidate,info
