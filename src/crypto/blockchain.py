import hashlib
import json
import datetime


class Block:
    def __init__(self, index, timestamp, prev_hash, transaction) -> None:
        self.index = index
        self.timestamp = timestamp
        self.prev_hash = prev_hash
        self.transaction = transaction
        self.diff = 4  # 難易度 diffを追加
        self.now_hash = self.calc_hash()
        self.nonce = None  # 採掘時に計算する対象nonceを追加

    def calc_hash(self) -> int:
        # ここではnow_hash(自身のブロック)を除いたデータをjsonに変換し
        jsoned_data = {
            'index': self.index,
            'timestamp': self.timestamp,
            'prev_hash': self.prev_hash,
            'transaction': self.transaction,
            'diff': self.diff
        }
        json_text = json.dumps(jsoned_data, sort_keys=True)
        return hashlib.sha256(json_text.encode('ascii')).hexdigest()

    def check(self, nonce):
        nonce_joined = self.now_hash+str(nonce)
        calced = hashlib.sha256(nonce_joined.encode('ascii')).hexdigest()
        if calced[:self.diff].count('0') == self.diff:
            return True
        else:
            return False

    def mining(self, append_transaction):
        nonce = 0
        self.transaction.append(append_transaction)
        # 報酬の好きな取引を一つ入れたあ後にハッシュ値を再計算、このハッシュ値にnonceを足した上にdiff桁まで0が続くものを探していきます。
        self.now_hash = self.calc_hash()
        while True:
            nonce_joined = self.now_hash+str(nonce)
            calced = hashlib.sha256(nonce_joined.encode('ascii')).hexdigest()
            if calced[:self.diff:].count('0') == self.diff:
                break
            nonce += 1
        return nonce


block_chain = []

block = Block(0, str(datetime.datetime.now()), '-', [])
append_transaction = {'paiza': 'genesis_block'}
nonce = block.mining(append_transaction)  # 最初のブロックの採掘
block.nonce = nonce
block_chain.append(block)

for i in range(5):
    block = Block(i+1, str(datetime.datetime.now()),
                  block_chain[i].now_hash, ["取引データ"])
    append_transaction = {'paiza': '採掘報酬ゲット' + str(i)}
    nonce = block.mining(append_transaction)
    block.nonce = nonce
    block_chain.append(block)

for block in block_chain:
    # nonceとブロックのハッシュ値を使ってブロックがルールを満たしているか検証
    nonce_joined = block.now_hash+str(block.nonce)
    calced = hashlib.sha256(nonce_joined.encode('ascii')).hexdigest()

    print("index =", block.index, "sha256(",
          block.now_hash, "+", block.nonce, ") =", calced)
    print('transaction', block.transaction)
