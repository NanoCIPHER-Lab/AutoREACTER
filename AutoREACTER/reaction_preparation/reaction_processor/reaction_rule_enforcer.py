"""
reaction_rule_enforcer.py

Utilities for enforcing reaction-rule availability during species-pool
propagation.

This module controls which reaction rules are active during each propagation
iteration. Base reaction rules are enabled at the start, while conditional
reaction rules can be enabled only after their required trigger reactions have
been successfully processed.

The rule enforcer does not detect functional groups, generate reaction
instances, or apply reactions. It only decides whether a reaction rule is
currently allowed and updates the active rule set based on completed reaction
history.
"""
from dataclasses import dataclass, field
from typing import Iterable
import json

@dataclass
class ReactionRuleEnforcer:
    base_rules: set[str] = field(default_factory=set)
    conditional_rules: list[dict] = field(default_factory=list)
    active_rules: set[str] = field(default_factory=set)
    completed_rules: set[str] = field(default_factory=set)
    
    def __init__(self):
        self.rule_json_path = "reaction_rules.json"
        self.load_reaction_rules()

    def load_reaction_rules(self):
        """Load reaction rules from a JSON file."""
        with open(self.rule_json_path, "r") as f:
            self.reaction_rules = json.load(f)