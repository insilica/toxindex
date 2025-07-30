#!/usr/bin/env python3
"""
Test script for real use case: "Is PFAO carcinogenic?"
"""
import sys
import os
import json
from pathlib import Path
from datetime import datetime

# Add the webserver directory to the Python path
sys.path.insert(0, str(Path(__file__).parent / "webserver"))

# Load environment variables
try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

def test_real_case():
    """Test the agent with a real use case"""
    print("🧪 Testing real use case: PFAO carcinogenicity\n")
    
    try:
        from webserver.tools.deeptox_agent import deeptox_agent
        from webserver.tools.toxicity_models import ChemicalToxicityAssessment
        
        # Real test query
        test_query = "Is PFAO carcinogenic?"
        
        print(f"Query: {test_query}")
        print("Running agent...")
        
        # Run the agent
        response = deeptox_agent.run(test_query)
        
        if isinstance(response.content, ChemicalToxicityAssessment):
            print("✅ Agent returned ChemicalToxicityAssessment object")
            
            # Convert to JSON
            json_output = response.content.model_dump()
            
            # Save to file
            timestamp = datetime.now().isoformat().replace(':', '-')
            filename = f"pfao_carcinogenicity_{timestamp}.json"
            
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(json_output, f, indent=2, ensure_ascii=False)
            
            print(f"✅ Saved JSON output to {filename}")
            
            # Analyze the results
            analyze_results(response.content)
            
            return True
            
        elif isinstance(response.content, dict):
            print("✅ Agent returned dict object")
            
            # Save to file
            timestamp = datetime.now().isoformat().replace(':', '-')
            filename = f"pfao_carcinogenicity_{timestamp}.json"
            
            with open(filename, 'w', encoding='utf-8') as f:
                json.dump(response.content, f, indent=2, ensure_ascii=False)
            
            print(f"✅ Saved JSON output to {filename}")
            
            # Try to convert to ChemicalToxicityAssessment for analysis
            try:
                assessment = ChemicalToxicityAssessment(**response.content)
                analyze_results(assessment)
            except Exception as e:
                print(f"⚠️  Could not convert to ChemicalToxicityAssessment: {e}")
                print("📋 Raw response structure:")
                print(json.dumps(response.content, indent=2))
            
            return True
            
        else:
            print(f"⚠️  Response type: {type(response.content)}")
            print(f"Content preview: {str(response.content)[:200]}...")
            
            # Try to convert to ChemicalToxicityAssessment if it's already one
            if hasattr(response.content, 'model_dump'):
                try:
                    json_output = response.content.model_dump()
                    
                    # Save to file
                    timestamp = datetime.now().isoformat().replace(':', '-')
                    filename = f"pfao_carcinogenicity_{timestamp}.json"
                    
                    with open(filename, 'w', encoding='utf-8') as f:
                        json.dump(json_output, f, indent=2, ensure_ascii=False)
                    
                    print(f"✅ Saved JSON output to {filename}")
                    
                    # Analyze the results
                    analyze_results(response.content)
                    
                    return True
                except Exception as e:
                    print(f"❌ Error processing response: {e}")
                    return False
            else:
                print(f"❌ Cannot process response type: {type(response.content)}")
                return False
            
    except Exception as e:
        print(f"❌ Error testing agent: {e}")
        import traceback
        traceback.print_exc()
        return False

def analyze_results(assessment):
    """Analyze the results of the assessment"""
    print("\n📊 Analyzing results...")
    
    # Check if we have beta distribution data
    try:
        beta_summary = assessment.get_beta_distribution_summary()
        print(f"📈 Beta Distribution Analysis:")
        print(f"   - Carcinogenicity probability: {beta_summary['probability']:.3f} ({beta_summary['probability']*100:.1f}%)")
        print(f"   - Confidence interval: [{beta_summary['confidence_interval']['lower']:.3f}, {beta_summary['confidence_interval']['upper']:.3f}]")
        print(f"   - Confidence level: {beta_summary['confidence_interval']['level']*100:.0f}%")
        print(f"   - Studies used: {beta_summary['studies_count']}")
        print(f"   - Overall confidence: {beta_summary['confidence']}")
        
        # Interpret the probability
        prob = beta_summary['probability']
        if prob < 0.1:
            risk_level = "Low"
        elif prob < 0.3:
            risk_level = "Moderate"
        elif prob < 0.5:
            risk_level = "High"
        else:
            risk_level = "Very High"
        
        print(f"   - Risk assessment: {risk_level} carcinogenicity risk")
        
    except Exception as e:
        print(f"❌ Error analyzing beta distribution: {e}")
    
    # Check reference quality
    try:
        quality_report = assessment.get_reference_quality_report()
        print(f"\n🔍 Reference Quality:")
        print(f"   - Total references: {quality_report['total_references']}")
        print(f"   - Valid references: {quality_report['valid_references']}")
        print(f"   - Invalid references: {quality_report['invalid_references']}")
        print(f"   - Unreachable references: {quality_report['unreachable_references']}")
        print(f"   - Validity rate: {quality_report['validity_rate']:.1%}")
        
        if quality_report['validity_rate'] < 0.5:
            print("   ⚠️  Low reference validity - many fake or inaccessible URLs")
        elif quality_report['validity_rate'] < 0.8:
            print("   ⚠️  Moderate reference validity - some issues with URLs")
        else:
            print("   ✅ Good reference validity")
            
    except Exception as e:
        print(f"❌ Error analyzing reference quality: {e}")
    
    # Check data completeness
    try:
        completeness = assessment.metadata.data_completeness
        print(f"\n📋 Data Completeness:")
        print(f"   - Overall score: {completeness.overall_score:.2f} ({completeness.overall_score*100:.0f}%)")
        print(f"   - Confidence level: {completeness.confidence_level}")
        
        if completeness.overall_score < 0.5:
            print("   ⚠️  Low data completeness - limited information available")
        elif completeness.overall_score < 0.8:
            print("   ⚠️  Moderate data completeness - some gaps in information")
        else:
            print("   ✅ High data completeness - comprehensive information")
            
    except Exception as e:
        print(f"❌ Error analyzing data completeness: {e}")
    
    # Check if we have chemical properties
    try:
        chem_props = assessment.chemical_properties
        print(f"\n🧪 Chemical Properties:")
        print(f"   - Description: {chem_props.description}")
        print(f"   - Properties analyzed: {len(chem_props.properties)}")
        print(f"   - Evidence sources: {len(chem_props.evidence)}")
        
    except Exception as e:
        print(f"❌ Error analyzing chemical properties: {e}")
    
    # Check if we have toxicity mechanisms
    try:
        tox_mech = assessment.toxicity_mechanisms
        print(f"\n⚗️  Toxicity Mechanisms:")
        print(f"   - Description: {tox_mech.description}")
        print(f"   - Mechanisms identified: {len(tox_mech.mechanisms)}")
        
    except Exception as e:
        print(f"❌ Error analyzing toxicity mechanisms: {e}")
    
    # Check if we have clinical evidence
    try:
        clin_evid = assessment.clinical_evidence
        print(f"\n🏥 Clinical Evidence:")
        print(f"   - Description: {clin_evid.description}")
        print(f"   - Studies analyzed: {len(clin_evid.studies)}")
        print(f"   - Treatment protocols: {len(clin_evid.treatment_protocols)}")
        
    except Exception as e:
        print(f"❌ Error analyzing clinical evidence: {e}")
    
    # Check risk factors
    try:
        risk_factors = assessment.risk_factors
        print(f"\n⚠️  Risk Factors:")
        print(f"   - High-risk groups: {len(risk_factors.high_risk_groups)}")
        print(f"   - Modifying factors: {len(risk_factors.modifying_factors)}")
        print(f"   - Preventive measures: {len(risk_factors.preventive_measures)}")
        
    except Exception as e:
        print(f"❌ Error analyzing risk factors: {e}")

def main():
    """Run the test"""
    print("🚀 Testing real use case: PFAO carcinogenicity\n")
    
    success = test_real_case()
    
    if success:
        print("\n🎉 Real use case test completed!")
        print("📁 Check the generated JSON file for the complete analysis")
        print("📊 Review the analysis above for key findings")
    else:
        print("\n❌ Test failed!")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main()) 