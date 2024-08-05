import threading
import time

class StoppableThread(threading.Thread):
    def __init__(self, target=None):
        super().__init__()
        self._stop_event = threading.Event()
        self._target = target

    def run(self):
        while not self._stop_event.is_set():
            if self._target:
                self._target()
            time.sleep(10)

    def stop(self):
        self._stop_event.set()